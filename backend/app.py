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
    
    try:
        # 从PubMed获取文章
        articles = get_pubmed_articles(query, time_period)
        
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

# 静态文件路由
@app.route('/<path:path>')
def static_files(path):
    return send_from_directory(app.static_folder, path)

if __name__ == '__main__':
    # 获取环境变量中的端口，如果没有则使用默认值5000
    port = int(os.environ.get('PORT', 5000))
    # 启动应用，在开发模式下启用调试
    app.run(host='0.0.0.0', port=port, debug=True)