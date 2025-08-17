import os
import time
import json
import logging
import requests
from dotenv import load_dotenv

# 配置日志
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger('OllamaAPI')

# 加载环境变量
load_dotenv()

class OllamaAPI:
    def __init__(self):
        # Ollama 服务配置 - 强制使用localhost以解决Windows连接问题
        # 忽略环境变量中的OLLAMA_HOST，直接使用localhost
        self.ollama_host = "localhost"  # 强制设置为localhost，不读取环境变量
        self.ollama_port = os.getenv("OLLAMA_PORT", "11434")
        self.model_name = "deepseek-r1:1.5b"
        self.model_path = os.getenv("OLLAMA_MODEL_PATH", "F:\\Ollama\\models")
        
        # 构建正确的API URL
        self.api_url = f"http://{self.ollama_host}:{self.ollama_port}/api/chat"
        
        # 验证Ollama服务状态
        self._check_ollama_service()

    def _check_ollama_service(self):
        """检查Ollama服务是否可用"""
        try:
            # 直接连接Ollama服务，不通过代理
            response = requests.get(f"http://{self.ollama_host}:{self.ollama_port}", timeout=10)
            if response.status_code == 200:
                logger.info("成功连接到Ollama服务")
            else:
                logger.warning(f"Ollama服务响应异常，状态码: {response.status_code}")
        except requests.exceptions.ConnectionError:
            logger.error("无法连接到Ollama服务，请确保服务已启动")
            logger.error(f"请检查Ollama服务是否在 {self.ollama_host}:{self.ollama_port} 运行")
        except Exception as e:
            logger.error(f"检查Ollama服务时发生错误: {str(e)}")

    def chat(self, messages, temperature=0.5, max_tokens=4096):
        """调用 Ollama 模型进行对话"""
        # 构建请求数据
        payload = {
            "model": self.model_name,
            "messages": messages,
            "stream": False,  # 非流式返回
            "options": {
                "temperature": temperature,
                "max_tokens": max_tokens,
                "num_ctx": 4096  # 明确设置上下文长度，适应低显存模式
            }
        }
        
        logger.debug(f"请求数据: {json.dumps(payload, ensure_ascii=False)[:200]}...")
        
        try:
            # 发送请求，不使用系统代理
            session = requests.Session()
            session.trust_env = False  # 禁用环境变量中的代理设置
            
            response = session.post(
                self.api_url,
                headers={"Content-Type": "application/json"},
                json=payload,
                timeout=300  # 进一步延长超时时间，适应低显存模式下的模型加载
            )
            
            if response.status_code != 200:
                raise Exception(f"API 调用失败: {response.text} (状态码: {response.status_code})")
            
            result = response.json()
            logger.debug(f"收到响应: {json.dumps(result, ensure_ascii=False)[:200]}...")
            
            # 提取回复内容
            if "message" in result and "content" in result["message"]:
                return result["message"]["content"]
            else:
                raise Exception("无法从响应中提取内容")
                
        except requests.exceptions.ConnectionError:
            error_msg = f"无法连接到Ollama服务，请确保服务已在 {self.ollama_host}:{self.ollama_port} 启动"
            logger.error(error_msg)
            raise Exception(error_msg)
        except Exception as e:
            logger.error(f"调用 Ollama API 失败: {str(e)}", exc_info=True)
            raise Exception(f"调用失败: {str(e)}")
    