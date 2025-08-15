// 主要JavaScript文件

// DOM元素
const searchInput = document.getElementById('search-input');
const searchButton = document.getElementById('search-button');
const timePeriod = document.getElementById('time-period');
const journalRank = document.getElementById('journal-rank');
const sortBy = document.getElementById('sort-by');
const resultsContainer = document.getElementById('results-container');
const resultCount = document.getElementById('result-count');
const loadingSpinner = document.getElementById('loading');
const articleTemplate = document.getElementById('article-template');

// API URL
const API_BASE_URL = '/api';

// 初始化
document.addEventListener('DOMContentLoaded', function() {
    // 绑定搜索按钮点击事件
    if (searchButton) {
        searchButton.addEventListener('click', performSearch);
    }
    
    // 绑定回车键搜索
    if (searchInput) {
        searchInput.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') {
                performSearch();
            }
        });
    }
    
    // 绑定排序变化事件
    if (sortBy) {
        sortBy.addEventListener('change', () => {
            if (resultsContainer.querySelectorAll('.article-card').length > 0) {
                performSearch();
            }
        });
    }
});

// 执行搜索
async function performSearch() {
    const query = searchInput.value.trim();
    
    if (!query) {
        alert('请输入搜索关键词');
        return;
    }
    
    // 显示加载动画
    showLoading(true);
    
    // 清空结果容器
    clearResults();
    
    try {
        // 构建查询参数
        const params = new URLSearchParams({
            query: query,
            time_period: timePeriod.value,
            sort: sortBy.value
        });
        
        // 发送API请求
        const response = await fetch(`${API_BASE_URL}/articles?${params}`);
        const data = await response.json();
        
        if (data.status === 'success') {
            // 显示结果数量
            resultCount.textContent = `(${data.data.length})`;
            
            // 根据期刊分区筛选
            let filteredArticles = data.data;
            if (journalRank.value !== 'all') {
                filteredArticles = filterArticlesByRank(data.data, journalRank.value);
                resultCount.textContent = `(${filteredArticles.length})`;
            }
            
            // 显示结果
            if (filteredArticles.length > 0) {
                displayArticles(filteredArticles);
            } else {
                showEmptyState('未找到符合条件的文献');
            }
        } else {
            showEmptyState('搜索出错: ' + data.message);
        }
    } catch (error) {
        console.error('搜索出错:', error);
        showEmptyState('搜索请求失败，请稍后重试');
    } finally {
        // 隐藏加载动画
        showLoading(false);
    }
}

// 根据期刊分区筛选文章
function filterArticlesByRank(articles, rankFilter) {
    return articles.filter(article => {
        if (!article.journal_metrics) return false;
        
        const { cas_division } = article.journal_metrics;
        
        switch (rankFilter) {
            case 'cas-1':
                return cas_division === '1区';
            case 'cas-2':
                return cas_division === '2区';
            default:
                return true;
        }
    });
}

// 显示文章列表
function displayArticles(articles) {
    articles.forEach(article => {
        const articleElement = createArticleElement(article);
        resultsContainer.appendChild(articleElement);
    });
}

// 创建文章元素
function createArticleElement(article) {
    const template = articleTemplate.content.cloneNode(true);
    const articleCard = template.querySelector('.article-card');
    
    // 设置标题
    template.querySelector('.article-title').textContent = article.title;
    
    // 设置作者
    if (article.authors && article.authors.length > 0) {
        template.querySelector('.authors').textContent = article.authors.join(', ');
    } else {
        template.querySelector('.authors').textContent = '作者信息不可用';
    }
    
    // 设置期刊信息
    if (article.journal) {
        const journalName = template.querySelector('.journal-name');
        const publishDate = template.querySelector('.publish-date');
        
        journalName.textContent = article.journal.title || '期刊未知';
        publishDate.textContent = article.journal.pub_date ? `(${article.journal.pub_date})` : '';
    }
    
    // 设置期刊指标
    if (article.journal_metrics) {
        const metrics = article.journal_metrics;
        
        // 隐藏JCR分区
        template.querySelector('.jcr-quartile').classList.add('hidden');
        
        if (metrics.cas_division) {
            const casElement = template.querySelector('.cas-division');
            casElement.textContent = `中科院: ${metrics.cas_division}`;
            casElement.classList.remove('hidden');
        } else {
            template.querySelector('.cas-division').classList.add('hidden');
        }
        
        // 隐藏影响因子
        template.querySelector('.impact-factor').classList.add('hidden');
    }
    
    // 设置摘要
    const abstractElement = template.querySelector('.article-abstract');
    if (article.abstract) {
        // 限制摘要长度
        const maxLength = 250;
        if (article.abstract.length > maxLength) {
            abstractElement.textContent = article.abstract.substring(0, maxLength) + '...';
        } else {
            abstractElement.textContent = article.abstract;
        }
    } else {
        abstractElement.textContent = '摘要不可用';
    }
    
    // 文章卡片不再需要点击事件
    // 移除了模态框，不再需要点击显示详情
    
    // 设置链接
    const doiLink = template.querySelector('.doi-link');
    if (article.doi) {
        doiLink.href = `https://doi.org/${article.doi}`;
        doiLink.target = "_blank"; // 确保在新标签页打开
        doiLink.onclick = (e) => {
            e.stopPropagation(); // 阻止事件冒泡，防止触发卡片点击
        };
    } else {
        doiLink.classList.add('hidden');
    }
    
    const pubmedLink = template.querySelector('.pubmed-link');
    if (article.pmid) {
        pubmedLink.href = `https://pubmed.ncbi.nlm.nih.gov/${article.pmid}/`;
        pubmedLink.target = "_blank"; // 确保在新标签页打开
        pubmedLink.onclick = (e) => {
            e.stopPropagation(); // 阻止事件冒泡，防止触发卡片点击
        };
    } else {
        pubmedLink.classList.add('hidden');
    }
    
    // 确保所有按钮点击不会触发卡片点击事件
    const allButtons = template.querySelectorAll('.article-actions a');
    allButtons.forEach(button => {
        button.addEventListener('click', (e) => {
            e.stopPropagation();
        });
    });
    
    // 不再需要移除查看详情按钮，因为模板中已经不存在该按钮
    
    return template;
}

// 模态框相关功能已移除

// 清空结果容器
function clearResults() {
    resultsContainer.innerHTML = '';
}

// 显示空状态
function showEmptyState(message) {
    const emptyState = document.createElement('div');
    emptyState.className = 'empty-state';
    emptyState.innerHTML = `
        <i class="fas fa-search fa-3x"></i>
        <p>${message}</p>
    `;
    resultsContainer.appendChild(emptyState);
}

// 显示/隐藏加载动画
function showLoading(show) {
    if (show) {
        loadingSpinner.classList.remove('hidden');
    } else {
        loadingSpinner.classList.add('hidden');
    }
}