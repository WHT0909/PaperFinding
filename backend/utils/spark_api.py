import os
import time
import json
import hashlib
import base64
import hmac
import logging
from dotenv import load_dotenv

# 配置日志以方便调试
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger('SparkAPI')

# 加载环境变量
load_dotenv()

class Spark4UltraAPI:
    def __init__(self):
        # 从环境变量获取配置
        self.appid = os.getenv("SPARK_APPID")
        self.api_secret = os.getenv("SPARK_API_SECRET")
        self.api_key = os.getenv("SPARK_API_KEY")
        
        # Spark4.0 Ultra 正确配置
        self.domain = "4.0Ultra"  # 4.0 Ultra专用domain
        self.api_version = "v4.0"  # API版本
        self.url = f"wss://spark-api.xf-yun.com/{self.api_version}/chat"
        
        # 验证配置完整性
        self._validate_config()

    def _validate_config(self):
        """验证配置是否完整"""
        missing = []
        if not self.appid:
            missing.append("SPARK_APPID")
        if not self.api_secret:
            missing.append("SPARK_API_SECRET")
        if not self.api_key:
            missing.append("SPARK_API_KEY")
            
        if missing:
            raise ValueError(f"缺少必要配置: {', '.join(missing)}")
        
        # 日志输出部分隐藏的密钥信息，方便调试但不泄露完整信息
        logger.debug(f"配置检查: appid={self.appid[:6]}..., api_key={self.api_key[:6]}...")

    def _create_signature(self, date_str):
        """生成签名 - 严格按照讯飞API文档实现"""
        # 签名基础字符串必须与请求路径完全匹配
        signature_origin = f"host: spark-api.xf-yun.com\ndate: {date_str}\nGET /{self.api_version}/chat HTTP/1.1"
        logger.debug(f"签名原始字符串: {signature_origin}")
        
        # 使用HMAC-SHA256计算签名
        signature_sha = hmac.new(
            self.api_secret.encode('utf-8'),
            signature_origin.encode('utf-8'),
            digestmod=hashlib.sha256
        ).digest()
        
        # Base64编码
        signature = base64.b64encode(signature_sha).decode(encoding='utf-8')
        logger.debug(f"生成的签名: {signature[:10]}...")  # 只显示部分签名
        return signature

    def _get_auth_headers(self):
        """获取认证头信息 - 确保日期格式正确"""
        # 生成严格符合RFC 1123格式的GMT时间
        date = time.strftime("%a, %d %b %Y %H:%M:%S GMT", time.gmtime())
        logger.debug(f"生成的日期: {date}")
        
        signature = self._create_signature(date)
        
        # 构建Authorization原始字符串
        authorization_origin = (
            f'api_key="{self.api_key}", '
            f'signature="{signature}", '
            f'date="{date}", '
            f'host="spark-api.xf-yun.com"'
        )
        
        # Base64编码Authorization
        authorization = base64.b64encode(authorization_origin.encode('utf-8')).decode(encoding='utf-8')
        logger.debug(f"生成的Authorization: {authorization[:10]}...")  # 只显示部分
        
        return {
            "authorization": authorization,
            "date": date,
            "host": "spark-api.xf-yun.com"
        }

    def chat(self, messages, temperature=0.5, max_tokens=4096, top_k=4):
        """调用Spark4.0 Ultra API进行对话"""
        import websockets
        import asyncio
        
        # 构建请求数据 - 严格匹配4.0 Ultra的要求
        data = {
            "header": {
                "app_id": self.appid,
                "uid": f"user_{int(time.time())}"  # 生成唯一用户ID
            },
            "parameter": {
                "chat": {
                    "domain": self.domain,
                    "temperature": temperature,
                    "max_tokens": max_tokens,
                    "top_k": top_k,  # 4.0 Ultra支持的参数
                    "response_format": "text"  # 明确指定响应格式
                }
            },
            "payload": {
                "message": {
                    "text": messages
                }
            }
        }
        
        logger.debug(f"请求数据: {json.dumps(data, ensure_ascii=False)[:200]}...")
        
        try:
            # 生成带签名的URL
            headers = self._get_auth_headers()
            
            # 正确编码URL参数
            url = f"{self.url}?authorization={headers['authorization']}&date={headers['date']}&host={headers['host']}"
            logger.debug(f"连接URL: {url[:100]}...")  # 只显示部分URL
            
            # 异步处理WebSocket
            async def _send_request():
                try:
                    # 增加超时设置
                    conn = websockets.connect(
                        url,
                        timeout=30,
                        ping_interval=10,
                        ping_timeout=5
                    )
                    
                    async with conn as websocket:
                        logger.debug("WebSocket连接成功")
                        await websocket.send(json.dumps(data))
                        logger.debug("请求已发送")
                        
                        full_response = ""
                        while True:
                            response = await websocket.recv()
                            logger.debug(f"收到响应: {response[:100]}...")
                            
                            response_json = json.loads(response)
                            
                            # 检查API返回的错误码
                            header = response_json.get("header", {})
                            if header.get("code") != 0:
                                raise Exception(f"API错误: {header.get('message')} (错误码: {header.get('code')})")
                            
                            # 提取回复内容
                            payload = response_json.get("payload", {})
                            choices = payload.get("choices", {})
                            text = choices.get("text", [])
                            if text:
                                full_response += text[0].get("content", "")
                            
                            # 检查是否结束
                            if choices.get("status") == 2:
                                logger.debug("响应完成")
                                break
                        
                        return full_response
                        
                except websockets.exceptions.InvalidStatusCode as e:
                    # 详细输出认证失败信息
                    logger.error(f"WebSocket认证失败: 状态码 {e.status_code}")
                    if e.status_code == 401:
                        logger.error("可能的原因: 签名错误、API密钥错误、服务未开通、IP限制")
                    raise Exception(f"认证失败 (状态码: {e.status_code})，请检查签名生成和API配置")
                except Exception as e:
                    logger.error(f"WebSocket通信错误: {str(e)}", exc_info=True)
                    raise Exception(f"通信错误: {str(e)}")
            
            # 执行异步请求
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            try:
                result = loop.run_until_complete(_send_request())
                return result
            finally:
                loop.close()
                
        except Exception as e:
            logger.error(f"调用星火API失败: {str(e)}", exc_info=True)
            raise Exception(f"调用失败: {str(e)}")
