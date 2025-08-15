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
    
    // 设置摘要
    if (article.abstract) {
        articleAbstract.innerHTML = formatAbstract(article.abstract);
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
}

// 格式化摘要文本
function formatAbstract(abstract) {
    if (!abstract || abstract.trim() === '') {
        return '<p class="no-abstract">摘要不可用</p>';
    }
    // 直接放进一个 <p> 标签，保持 PubMed 返回的原始文本
    return `<p>${abstract.trim()}</p>`;
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