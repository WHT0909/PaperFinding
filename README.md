# PaperFinding - 文献检索软件

## 项目简介

PaperFinding是一个帮助研究人员查看热门领域（如AI、脑机接口BCI等）最新文献的检索软件。该软件通过集成多个学术数据库的API，提供友好的用户界面，帮助研究人员快速获取高质量的学术文献。

## 主要功能

1. 从多个数据库（如PubMed等）获取最新文献
2. 提供友好的用户界面，便于交互
3. 支持按发表时间筛选文献（如一年内、一个月内等）
4. 显示文章的中科院分区、JCR分区等质量指标
5. 提供文章的标题、作者、摘要等详细信息

## 技术栈

- 前端：HTML + CSS + JavaScript
- 后端：Python + Flask
- 数据来源：PubMed API等学术数据库

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

我的本地 deepseek 模型文件在路径 "F:\Ollama\models" 下。以管理员方式打开终端，使用`ollama server`指令启动 ollama；打开浏览器访问 http://localhost:11434，确认能看到 "Ollama is running" 后，进入下一步骤

注：如果有报错提示端口冲突 / 被占用，请先手动退出 ollama 应用程序

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
├── frontend/           # 前端文件
│   ├── index.html      # 主页面
│   ├── css/            # 样式文件
│   └── js/             # JavaScript文件
├── backend/            # 后端文件
│   ├── api/            # API接口
│   ├── database/       # 数据库操作
│   └── utils/          # 工具函数
├── venv/               # 虚拟环境
└── requirements.txt    # 依赖列表
```

## API密钥配置

部分数据库API可能需要密钥才能访问。请在项目根目录创建`.env`文件，并按以下格式添加您的API密钥：

```
PUBMED_API_KEY=your_pubmed_api_key
# 其他API密钥
```

## 贡献指南

欢迎贡献代码或提出建议！请通过Issue或Pull Request参与项目开发。

## 许可证

[MIT License](LICENSE)