import re
from datetime import datetime, timedelta

def format_date(date_str):
    """
    将各种格式的日期字符串转换为标准格式 (YYYY-MM-DD)
    
    参数:
    - date_str: 日期字符串
    
    返回:
    - 标准格式的日期字符串，如果无法解析则返回原字符串
    """
    try:
        # 尝试解析常见的日期格式
        formats = [
            '%Y-%m-%d', '%Y/%m/%d', '%d-%m-%Y', '%d/%m/%Y',
            '%b %d, %Y', '%B %d, %Y', '%Y-%m', '%Y/%m',
            '%Y'
        ]
        
        for fmt in formats:
            try:
                date_obj = datetime.strptime(date_str, fmt)
                return date_obj.strftime('%Y-%m-%d')
            except ValueError:
                continue
        
        # 如果上述格式都不匹配，尝试使用正则表达式提取年月日
        year_pattern = r'\b(19|20)\d{2}\b'
        month_pattern = r'\b(0?[1-9]|1[0-2])\b'
        day_pattern = r'\b(0?[1-9]|[12]\d|3[01])\b'
        
        year_match = re.search(year_pattern, date_str)
        month_match = re.search(month_pattern, date_str)
        day_match = re.search(day_pattern, date_str)
        
        if year_match:
            year = year_match.group()
            month = month_match.group() if month_match else '01'
            day = day_match.group() if day_match else '01'
            
            # 确保月和日是两位数
            month = month.zfill(2)
            day = day.zfill(2)
            
            return f"{year}-{month}-{day}"
        
        # 如果无法解析，返回原字符串
        return date_str
    
    except Exception:
        return date_str

def calculate_date_range(time_period):
    """
    根据时间段计算日期范围
    
    参数:
    - time_period: 时间段字符串，例如 '1year', '6months', '30days'
    
    返回:
    - 开始日期和结束日期的元组 (start_date, end_date)，格式为 YYYY-MM-DD
    """
    today = datetime.now()
    end_date = today.strftime('%Y-%m-%d')
    
    # 解析时间段
    if time_period.endswith('year') or time_period.endswith('years'):
        years = int(time_period.split('year')[0])
        start_date = today - timedelta(days=365 * years)
    elif time_period.endswith('month') or time_period.endswith('months'):
        months = int(time_period.split('month')[0])
        start_date = today - timedelta(days=30 * months)
    elif time_period.endswith('day') or time_period.endswith('days'):
        days = int(time_period.split('day')[0])
        start_date = today - timedelta(days=days)
    else:
        # 默认为1年
        start_date = today - timedelta(days=365)
    
    return start_date.strftime('%Y-%m-%d'), end_date

def truncate_text(text, max_length=200):
    """
    截断文本到指定长度，并添加省略号
    
    参数:
    - text: 要截断的文本
    - max_length: 最大长度
    
    返回:
    - 截断后的文本
    """
    if not text or len(text) <= max_length:
        return text
    
    return text[:max_length] + '...'

def clean_html(html_text):
    """
    清除HTML标签
    
    参数:
    - html_text: 包含HTML标签的文本
    
    返回:
    - 清除HTML标签后的纯文本
    """
    if not html_text:
        return ''
    
    # 使用正则表达式移除HTML标签
    clean_text = re.sub(r'<.*?>', '', html_text)
    # 替换多个空格为单个空格
    clean_text = re.sub(r'\s+', ' ', clean_text)
    # 去除首尾空格
    clean_text = clean_text.strip()
    
    return clean_text