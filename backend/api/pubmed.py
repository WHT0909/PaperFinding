import os
import requests
from datetime import datetime, timedelta
from Bio import Entrez
from dotenv import load_dotenv


# 获取API密钥
load_dotenv()          # 会把 .env 读进环境变量
API_KEY = os.getenv("PUBMED_API_KEY", "")
Entrez.email = os.getenv("Entrez.email", "")
def format_time_period(time_period):
    """
    将时间段转换为PubMed API可接受的格式
    例如：'1year', '6months', '30days'
    """
    today = datetime.now()
    
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
    
    # 格式化为YYYY/MM/DD格式
    return start_date.strftime("%Y/%m/%d")

def get_pubmed_articles(query, time_period='1year', max_results=100):
    """
    从PubMed获取文章
    
    参数:
    - query: 搜索查询
    - time_period: 时间段 (例如 '1year', '6months', '30days')
    - max_results: 最大结果数
    
    返回:
    - 文章列表
    """
    try:
        # 格式化时间段
        date_range = format_time_period(time_period)
        
        # 构建查询
        search_query = f"{query} AND (\"{date_range}\"[Date - Publication] : \"3000\"[Date - Publication])"
        
        # 使用Entrez API搜索PubMed
        if API_KEY:
            handle = Entrez.esearch(db="pubmed", term=search_query, retmax=max_results, sort="pub_date", api_key=API_KEY)
        else:
            handle = Entrez.esearch(db="pubmed", term=search_query, retmax=max_results, sort="pub_date")
            
        record = Entrez.read(handle)
        handle.close()
        
        # 获取文章ID列表
        id_list = record["IdList"]
        
        if not id_list:
            return []
        
        # 获取文章详情
        if API_KEY:
            handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml", api_key=API_KEY)
        else:
            handle = Entrez.efetch(db="pubmed", id=id_list, retmode="xml")
            
        records = Entrez.read(handle)
        handle.close()
        
        # 解析文章信息
        articles = []
        for record in records["PubmedArticle"]:
            article_data = {}
            article = record["MedlineCitation"]["Article"]
            
            # 标题
            article_data["title"] = article["ArticleTitle"]
            
            # 作者
            if "AuthorList" in article:
                authors = []
                for author in article["AuthorList"]:
                    if "LastName" in author and "ForeName" in author:
                        authors.append(f"{author['LastName']} {author['ForeName']}")
                    elif "CollectiveName" in author:
                        authors.append(author["CollectiveName"])
                article_data["authors"] = authors
            
            # 期刊
            if "Journal" in article:
                journal = article["Journal"]
                journal_info = {}
                if "Title" in journal:
                    journal_info["title"] = journal["Title"]
                if "ISSN" in journal:
                    journal_info["issn"] = journal["ISSN"]
                if "JournalIssue" in journal:
                    issue = journal["JournalIssue"]
                    if "Volume" in issue:
                        journal_info["volume"] = issue["Volume"]
                    if "Issue" in issue:
                        journal_info["issue"] = issue["Issue"]
                    if "PubDate" in issue:
                        pub_date = issue["PubDate"]
                        date_parts = []
                        for key in ["Year", "Month", "Day"]:
                            if key in pub_date:
                                date_parts.append(pub_date[key])
                        journal_info["pub_date"] = "-".join(date_parts)
                article_data["journal"] = journal_info
            
            # 摘要
            if "Abstract" in article and "AbstractText" in article["Abstract"]:
                abstract_parts = article["Abstract"]["AbstractText"]
                if isinstance(abstract_parts, list):
                    abstract = " ".join([str(part) for part in abstract_parts])
                else:
                    abstract = str(abstract_parts)
                article_data["abstract"] = abstract
            
            # DOI
            if "ELocationID" in article:
                for location in article["ELocationID"]:
                    if location.attributes["EIdType"] == "doi":
                        article_data["doi"] = str(location)
            
            # PubMed ID
            article_data["pmid"] = record["MedlineCitation"]["PMID"]
            
            # 添加发表日期用于排序（从PubDate中提取）
            if "journal" in article_data and "pub_date" in article_data["journal"]:
                article_data["sort_date"] = article_data["journal"]["pub_date"]
            else:
                article_data["sort_date"] = "1900-01-01"  # 默认日期
            
            articles.append(article_data)
        
        # 按发表日期从新到旧排序
        def parse_date(date_str):
            """解析日期字符串，返回可比较的日期对象"""
            try:
                # 处理不同的日期格式
                if '-' in date_str:
                    parts = date_str.split('-')
                    if len(parts) >= 1:
                        year = int(parts[0])
                        month = int(parts[1]) if len(parts) > 1 else 1
                        day = int(parts[2]) if len(parts) > 2 else 1
                        return datetime(year, month, day)
                else:
                    # 只有年份的情况
                    return datetime(int(date_str), 1, 1)
            except:
                return datetime(1900, 1, 1)
        
        # 按日期排序（从新到旧）
        articles.sort(key=lambda x: parse_date(x.get('sort_date', '1900-01-01')), reverse=True)
        
        return articles
    
    except Exception as e:
        print(f"Error fetching PubMed articles: {e}")
        return []