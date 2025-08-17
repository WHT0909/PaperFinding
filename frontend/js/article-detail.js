// 文章详情页面JavaScript

// DOM元素
const loadingElement = document.getElementById('loading');
const errorElement = document.getElementById('error-message');
const errorText = document.getElementById('error-text');
const articleDetail = document.getElementById('article-detail');
const articleTitle = document.getElementById('article-title');
const articleAuthors = document.getElementById('article-authors');
const journalName = document.getElementById('journal-name');
const publishDate = document.getElementById('publish-date');
const journalMetrics = document.getElementById('journal-metrics');
const casDivision = document.getElementById('cas-division');
const articleAbstract = document.getElementById('article-abstract');
const doiLink = document.getElementById('doi-link');
const pubmedLink = document.getElementById('pubmed-link');
// 添加对话相关DOM元素
const chatContainer = document.getElementById('chat-container');
const messagesContainer = document.getElementById('messages-container');
const userInput = document.getElementById('user-input');
const sendButton = document.getElementById('send-button');
let currentArticle = null;

// API URL
const API_BASE_URL = '/api';

// 初始化
document.addEventListener('DOMContentLoaded', function() {
    // 从URL参数获取PMID
    const urlParams = new URLSearchParams(window.location.search);
    const pmid = urlParams.get('pmid');
    
    if (!pmid) {
        showError('缺少文章ID参数');
        return;
    }
    
    // 加载文章详情
    loadArticleDetail(pmid);
});

// 加载文章详情
async function loadArticleDetail(pmid) {
    showLoading(true);
    hideError();
    hideArticleDetail();
    
    try {
        const response = await fetch(`${API_BASE_URL}/article/${pmid}`);
        const data = await response.json();
        
        if (data.status === 'success') {
            displayArticleDetail(data.data);
        } else {
            showError(data.message || '加载文章详情失败');
        }
    } catch (error) {
        console.error('加载文章详情出错:', error);
        showError('网络请求失败，请稍后重试');
    } finally {
        showLoading(false);
    }
}

// 显示文章详情
function displayArticleDetail(article) {
    // 设置标题
    articleTitle.textContent = article.title || '标题不可用';
    document.title = `${article.title || '文章详情'} - PaperFinding`;
    
    // 设置作者
    if (article.authors && article.authors.length > 0) {
        articleAuthors.textContent = article.authors.join(', ');
    } else {
        articleAuthors.textContent = '作者信息不可用';
    }
    
    // 设置期刊信息
    if (article.journal) {
        journalName.textContent = article.journal.title || '期刊未知';
        publishDate.textContent = article.journal.pub_date ? `发表日期: ${article.journal.pub_date}` : '';
    } else {
        journalName.textContent = '期刊信息不可用';
        publishDate.textContent = '';
    }
    
    // 设置期刊指标
    if (article.journal_metrics && article.journal_metrics.cas_division) {
        casDivision.textContent = `中科院分区: ${article.journal_metrics.cas_division}`;
        casDivision.classList.remove('hidden');
    } else {
        casDivision.classList.add('hidden');
    }
    
    // 设置摘要（支持Markdown）
    if (article.abstract) {
        articleAbstract.innerHTML = marked.parse(article.abstract);
    } else {
        articleAbstract.innerHTML = '<p class="no-abstract">摘要不可用</p>';
    }
    
    // 设置链接
    if (article.doi) {
        doiLink.href = `https://doi.org/${article.doi}`;
        doiLink.classList.remove('hidden');
    } else {
        doiLink.classList.add('hidden');
    }
    
    if (article.pmid) {
        pubmedLink.href = `https://pubmed.ncbi.nlm.nih.gov/${article.pmid}/`;
    }
    
    // 显示文章详情
    showArticleDetail();

    // 初始化对话功能
    initChat(article);
}

// 显示加载动画
function showLoading(show) {
    if (show) {
        loadingElement.classList.remove('hidden');
    } else {
        loadingElement.classList.add('hidden');
    }
}

// 显示错误信息
function showError(message) {
    errorText.textContent = message;
    errorElement.classList.remove('hidden');
}

// 隐藏错误信息
function hideError() {
    errorElement.classList.add('hidden');
}

// 显示文章详情
function showArticleDetail() {
    articleDetail.classList.remove('hidden');
}

// 隐藏文章详情
function hideArticleDetail() {
    articleDetail.classList.add('hidden');
}

// 初始化对话功能
function initChat(article) {
    currentArticle = article;
    bindChatEvents();
    
    // 发送初始请求，让AI分析文章
    setTimeout(() => {
        addMessage('bot', '正在分析文章内容，请稍候...');
        sendAnalysisRequest();
    }, 1000);
}

// 绑定对话事件
function bindChatEvents() {
    // 发送按钮点击事件
    sendButton.addEventListener('click', sendMessage);
    
    // 回车键发送消息
    userInput.addEventListener('keypress', (e) => {
        if (e.key === 'Enter' && !e.shiftKey) {
            e.preventDefault();
            sendMessage();
        }
    });
}

// 发送消息
async function sendMessage() {
    const message = userInput.value.trim();
    if (!message || !currentArticle) return;
    
    // 添加用户消息
    addMessage('user', message);
    userInput.value = '';
    
    // 显示正在输入状态
    const typingIndicator = addTypingIndicator();
    
    try {
        // 发送请求到后端
        const response = await fetch(`${API_BASE_URL}/article/${currentArticle.pmid}/chat`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                message: message,
                article: {
                    title: currentArticle.title,
                    abstract: currentArticle.abstract,
                    authors: currentArticle.authors,
                    journal: currentArticle.journal
                }
            })
        });
        
        const data = await response.json();
        
        // 移除正在输入状态
        removeTypingIndicator(typingIndicator);
        
        if (data.status === 'success') {
            // 假设后端返回思考过程和最终结果
            if (data.thought_process && data.thought_process.length > 0) {
                // 逐步显示思考过程
                for (const thought of data.thought_process) {
                    await typeMessage('thinking', thought, 50); // 50ms延迟，模拟实时思考
                    await new Promise(resolve => setTimeout(resolve, 300)); // 思考间隔
                }
            }
            // 显示最终回复（支持Markdown）
            addMarkdownMessage('bot', data.response);
        } else {
            addMessage('bot', `抱歉，处理请求时出错: ${data.message}`);
        }
    } catch (error) {
        console.error('对话请求失败:', error);
        removeTypingIndicator(typingIndicator);
        addMessage('bot', '抱歉，无法连接到服务，请稍后重试');
    }
}

// 发送文章分析请求
async function sendAnalysisRequest() {
    if (!currentArticle) return;
    
    try {
        const response = await fetch(`${API_BASE_URL}/article/${currentArticle.pmid}/analyze`, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                article: {
                    title: currentArticle.title,
                    abstract: currentArticle.abstract,
                    authors: currentArticle.authors,
                    journal: currentArticle.journal
                }
            })
        });
        
        const data = await response.json();
        
        // 清除"正在分析"消息
        messagesContainer.lastChild.remove();
        
        if (data.status === 'success') {
            // 假设后端返回思考过程和最终结果
            if (data.thought_process && data.thought_process.length > 0) {
                // 逐步显示思考过程
                for (const thought of data.thought_process) {
                    await typeMessage('thinking', thought, 50); // 50ms延迟，模拟实时思考
                    await new Promise(resolve => setTimeout(resolve, 300)); // 思考间隔
                }
            }
            // 显示最终分析结果（支持Markdown）
            addMarkdownMessage('bot', data.analysis);
            addMessage('bot', '我已分析这篇文章，您有什么问题想问我吗？');
        } else {
            addMessage('bot', `分析文章时出错: ${data.message}`);
        }
    } catch (error) {
        console.error('分析请求失败:', error);
        // 清除"正在分析"消息
        messagesContainer.lastChild.remove();
        addMessage('bot', '抱歉，分析文章时出现错误，请稍后重试');
    }
}

// 添加普通文本消息到对话容器
function addMessage(sender, text) {
    const messageDiv = document.createElement('div');
    messageDiv.className = `message ${sender}`;
    
    const contentDiv = document.createElement('div');
    contentDiv.className = 'message-content';
    contentDiv.textContent = text;
    
    messageDiv.appendChild(contentDiv);
    messagesContainer.appendChild(messageDiv);
    
    // 滚动到底部
    scrollToBottom();
}

// 添加Markdown格式消息到对话容器
function addMarkdownMessage(sender, markdownText) {
    const messageDiv = document.createElement('div');
    messageDiv.className = `message ${sender}`;
    
    const contentDiv = document.createElement('div');
    contentDiv.className = 'message-content markdown-content';
    contentDiv.innerHTML = marked.parse(markdownText); // 使用marked.js解析Markdown
    
    messageDiv.appendChild(contentDiv);
    messagesContainer.appendChild(messageDiv);
    
    // 滚动到底部
    scrollToBottom();
}

// 打字机效果显示消息（思考过程）
async function typeMessage(sender, text, delay = 50) {
    const messageDiv = document.createElement('div');
    messageDiv.className = `message ${sender}`;
    
    const contentDiv = document.createElement('div');
    contentDiv.className = 'message-content typewriter-text';
    contentDiv.textContent = '';
    
    messageDiv.appendChild(contentDiv);
    messagesContainer.appendChild(messageDiv);
    scrollToBottom();
    
    // 逐字显示文本
    for (let i = 0; i < text.length; i++) {
        contentDiv.textContent += text.charAt(i);
        scrollToBottom();
        await new Promise(resolve => setTimeout(resolve, delay));
    }
    
    // 移除打字机效果的光标
    contentDiv.classList.remove('typewriter-text');
    return messageDiv;
}

// 添加正在输入指示器
function addTypingIndicator() {
    const indicatorDiv = document.createElement('div');
    indicatorDiv.className = 'message bot typing-indicator';
    indicatorDiv.innerHTML = '<span></span><span></span><span></span>';
    
    messagesContainer.appendChild(indicatorDiv);
    scrollToBottom();
    
    return indicatorDiv;
}

// 移除正在输入指示器
function removeTypingIndicator(indicator) {
    if (indicator && indicator.parentNode === messagesContainer) {
        messagesContainer.removeChild(indicator);
    }
}

// 滚动到消息容器底部
function scrollToBottom() {
    messagesContainer.scrollTop = messagesContainer.scrollHeight;
}
