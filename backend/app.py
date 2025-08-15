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

# 主页路由 - 提供前端页面
@app.route('/')
def index():
    return send_from_directory(app.static_folder, 'index.html')

# API路由 - 获取文献
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

# API路由 - 获取单篇文章详情
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

# 静态文件路由
@app.route('/<path:path>')
def static_files(path):
    return send_from_directory(app.static_folder, path)

if __name__ == '__main__':
    # 获取环境变量中的端口，如果没有则使用默认值5000
    port = int(os.environ.get('PORT', 5000))
    # 启动应用，在开发模式下启用调试
    app.run(host='0.0.0.0', port=port, debug=True)