# PaperFinding - 文献检索系统

## 项目简介

PaperFinding是一个帮助研究人员查看生物医学热门领域（如生物医学工程BME、脑机接口BCI等）最新文献的文章检索系统。该软件通过集成PubMed学术数据库的API，提供友好的用户界面，帮助研究人员快速获取高质量的学术文献。

## 主要功能

1. 从PubMed获取最新文献
2. 提供友好的用户界面，便于交互
3. 支持按发表时间筛选文献（如一年内、一个月内等）
4. 显示文章的质量指标，如中科院分区（2025年最新分区）
5. 提供文章的标题、作者、摘要等详细信息

## 技术栈

- 前端：HTML + CSS + JavaScript
- 后端：Python + Flask
- 数据来源：PubMed API

## 安装与使用

### 环境要求

- Python 3.7+
- 虚拟环境（推荐）

### 安装步骤

1. 克隆仓库
```
git clone [仓库地址]
cd PaperFinding
```

2. 创建并激活虚拟环境
```
python -m venv venv
# Windows
venv\Scripts\activate
# Linux/Mac
source venv/bin/activate
```

3. 安装依赖
```
pip install -r requirements.txt
```

4. 启动 ollama

- 以管理员方式打开终端，使用`ollama server`指令启动 ollama
- 打开浏览器访问 http://localhost:11434
- 确认能看到 "Ollama is running" 后，进入下一步骤

注1：如果有报错提示端口冲突 / 被占用，请先手动退出 ollama 应用程序

注2：代码中默认的模型文件存储在路径"F:\Ollama\models"下

5. 运行应用
```
python backend/app.py
```

6. 在浏览器中访问
```
http://localhost:5000
```

## 项目结构

```
PaperFinding/
├── frontend/                    # 前端文件
│   ├── index.html               # 主页（搜索页面）
│   ├── article-detail.html      # 文章详情页
│   ├── css/                     # 样式文件
│   │   └── style.css            # 全局样式
│   │   
│   └── js/                      # JavaScript文件
│       ├── main.js              # 主页交互逻辑
│       └── article-detail.js    # 文章详情页交互逻辑
│         
├── backend/                     # 后端文件
│   ├── app.py                   # 应用入口
│   ├── api/                     # API接口
│   │   ├── __init__.py
│   │   ├── journal_info.py      # 文章相关接口
│   │   └── pubmed.py            # PubMed相关接口
│   ├── database/                # 数据库操作
│   │   ├── __init__.py
│   │   ├── cas_data.csv         # 2025年中科院分区数据
│   │   └── jcr_data.csv         # jcr分区数据（暂未使用）
│   └── utils/                   # 工具函数
│       ├── __init__.py
│       ├── ollama_api.py        # Ollama 接口
│       └── helpers.py           # 通用辅助函数
│
├── venv/                        # Python虚拟环境
├── requirements.txt             # Python依赖列表
├── .gitignore                   # Git忽略文件
└── README.md                    # 项目说明文档
```

## API密钥配置

部分数据库API可能需要密钥才能访问。请在项目根目录创建`.env`文件，并按以下格式添加您的API密钥：

```
# PubMed配置
PUBMED_API_KEY=your_pubmed_api_key
# Ollama 配置
OLLAMA_HOST=localhost
OLLAMA_PORT=11434
OLLAMA_MODEL_PATH=F:\Ollama\models（如果您的模型文件在其他位置，请修改此处）
```

## 贡献指南

欢迎贡献代码或提出建议！请通过Issue或Pull Request参与项目开发。

## 许可证

[MIT License](LICENSE)
