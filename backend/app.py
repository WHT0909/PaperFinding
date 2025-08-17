import os
from flask import Flask, jsonify, request, send_from_directory
from flask_cors import CORS
from dotenv import load_dotenv

# 加载环境变量
load_dotenv()

# 创建Flask应用
app = Flask(__name__, static_folder='../frontend')
CORS(app)  # 启用CORS以允许前端访问API

# 导入API模块
from api.pubmed import get_pubmed_articles
from api.journal_info import get_journal_metrics
# 导入Spark4.0 Ultra API工具类
from utils.spark_api import Spark4UltraAPI  # 修改导入的类名

# 初始化星火API客户端 - 使用Spark4.0 Ultra
spark_api = Spark4UltraAPI()  # 修改类名

# 主页路由 - 提供前端页面
@app.route('/')
def index():
    return send_from_directory(app.static_folder, 'index.html')

# API路由 - 获取文献（保持不变）
@app.route('/api/articles', methods=['GET'])
def get_articles():
    # 获取查询参数
    query = request.args.get('query', '')
    time_period = request.args.get('time_period', '1year')  # 默认为1年内
    max_results = request.args.get('max_results', 50)  # 默认最大50篇
    
    try:
        # 从PubMed获取文章
        articles = get_pubmed_articles(query, time_period, int(max_results))
        
        # 获取期刊指标信息
        for article in articles:
            if 'journal' in article:
                metrics = get_journal_metrics(article['journal'])
                article['journal_metrics'] = metrics
        
        return jsonify({
            'status': 'success',
            'data': articles
        })
    except Exception as e:
        return jsonify({
            'status': 'error',
            'message': str(e)
        }), 500

# API路由 - 获取单篇文章详情（保持不变）
@app.route('/api/article/<pmid>', methods=['GET'])
def get_article_detail(pmid):
    try:
        from Bio import Entrez
        import os
        
        API_KEY = os.getenv("PUBMED_API_KEY", "")
        
        # 获取文章详情
        if API_KEY:
            handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml", api_key=API_KEY)
        else:
            handle = Entrez.efetch(db="pubmed", id=pmid, retmode="xml")
            
        records = Entrez.read(handle)
        handle.close()
        
        if not records["PubmedArticle"]:
            return jsonify({
                'status': 'error',
                'message': '文章未找到'
            }), 404
        
        record = records["PubmedArticle"][0]
        article = record["MedlineCitation"]["Article"]
        
        # 解析详细信息
        article_detail = {
            'pmid': pmid,
            'title': article["ArticleTitle"],
            'authors': [],
            'journal': {},
            'abstract': '',
            'doi': ''
        }
        
        # 作者信息
        if "AuthorList" in article:
            for author in article["AuthorList"]:
                if "LastName" in author and "ForeName" in author:
                    article_detail["authors"].append(f"{author['LastName']} {author['ForeName']}")
                elif "CollectiveName" in author:
                    article_detail["authors"].append(author["CollectiveName"])
        
        # 期刊信息
        if "Journal" in article:
            journal = article["Journal"]
            if "Title" in journal:
                article_detail["journal"]["title"] = journal["Title"]
            if "JournalIssue" in journal:
                issue = journal["JournalIssue"]
                if "PubDate" in issue:
                    pub_date = issue["PubDate"]
                    date_parts = []
                    for key in ["Year", "Month", "Day"]:
                        if key in pub_date:
                            date_parts.append(pub_date[key])
                    article_detail["journal"]["pub_date"] = "-".join(date_parts)
        
        # 完整摘要
        if "Abstract" in article and "AbstractText" in article["Abstract"]:
            abstract_parts = article["Abstract"]["AbstractText"]
            if isinstance(abstract_parts, list):
                abstract = " ".join([str(part) for part in abstract_parts])
            else:
                abstract = str(abstract_parts)
            article_detail["abstract"] = abstract
        
        # DOI
        if "ELocationID" in article:
            for location in article["ELocationID"]:
                if location.attributes["EIdType"] == "doi":
                    article_detail["doi"] = str(location)
        
        # 获取期刊指标信息
        if article_detail['journal']:
            metrics = get_journal_metrics(article_detail['journal'])
            article_detail['journal_metrics'] = metrics
        
        return jsonify({
            'status': 'success',
            'data': article_detail
        })
        
    except Exception as e:
        return jsonify({
            'status': 'error',
            'message': str(e)
        }), 500

# API路由 - 分析文章（使用Spark4.0 Ultra的增强能力）
@app.route('/api/article/<pmid>/analyze', methods=['POST'])
def analyze_article(pmid):
    try:
        data = request.json
        article = data.get('article', {})
        
        if not article or not article.get('title'):
            return jsonify({
                'status': 'error',
                'message': '缺少文章信息'
            }), 400
        
        # 构建分析提示 - 利用Spark4.0 Ultra的增强分析能力
        prompt = f"""请深入分析以下学术文章并提供全面总结:
        
        标题: {article.get('title', '无标题')}
        作者: {', '.join(article.get('authors', ['未知']))}
        期刊: {article.get('journal', {}).get('title', '未知期刊')}
        发表日期: {article.get('journal', {}).get('pub_date', '未知日期')}
        
        摘要: {article.get('abstract', '无摘要')}
        
        请从以下几个方面进行分析:
        1. 研究背景与意义：该研究要解决什么问题？为什么这个问题重要？
        2. 核心方法：作者采用了什么研究方法或技术手段？
        3. 关键发现：研究得出了哪些重要结果或发现？
        4. 研究局限性：该研究存在哪些不足或局限性？
        5. 未来研究方向：基于该研究，未来可以开展哪些相关研究？
        
        请用专业、简洁的中文回答，总字数控制在500字左右。
        """
        
        # 调用Spark4.0 Ultra API，使用较低的temperature获得更稳定的分析结果
        messages = [
            {"role": "system", "content": "你是一位专业的学术文献分析专家，擅长深入解析各类学术论文并提取关键信息。"},
            {"role": "user", "content": prompt}
        ]
        
        analysis = spark_api.chat(messages, temperature=0.2, max_tokens=2048)
        
        return jsonify({
            'status': 'success',
            'analysis': analysis
        })
        
    except Exception as e:
        return jsonify({
            'status': 'error',
            'message': str(e)
        }), 500

# API路由 - 文章对话（保持不变，受益于Spark4.0 Ultra的增强能力）
@app.route('/api/article/<pmid>/chat', methods=['POST'])
def chat_about_article(pmid):
    try:
        data = request.json
        message = data.get('message', '')
        article = data.get('article', {})
        
        if not message:
            return jsonify({
                'status': 'error',
                'message': '请输入问题'
            }), 400
        
        # 构建对话提示
        prompt = f"""基于以下文章信息回答用户问题:
        
        标题: {article.get('title', '无标题')}
        作者: {', '.join(article.get('authors', ['未知']))}
        期刊: {article.get('journal', {}).get('title', '未知期刊')}
        
        摘要: {article.get('abstract', '无摘要')}
        
        用户问题: {message}
        """
        
        # 调用Spark4.0 Ultra API
        messages = [
            {"role": "system", "content": "你是一位专业的学术文献分析助手，根据提供的文章信息详细回答用户问题，用中文回答。"},
            {"role": "user", "content": prompt}
        ]
        
        answer = spark_api.chat(messages, temperature=0.5, max_tokens=4096)
        
        return jsonify({
            'status': 'success',
            'response': answer
        })
        
    except Exception as e:
        return jsonify({
            'status': 'error',
            'message': str(e)
        }), 500

# 静态文件路由
@app.route('/<path:path>')
def static_files(path):
    return send_from_directory(app.static_folder, path)

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=True)
