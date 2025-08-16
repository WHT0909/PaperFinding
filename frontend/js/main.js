// 主要JavaScript文件

// DOM元素
const searchInput = document.getElementById('search-input');
const searchButton = document.getElementById('search-button');
const timePeriod = document.getElementById('time-period');
const journalRank = document.getElementById('journal-rank');
const maxResults = document.getElementById('max-results');
const sortBy = document.getElementById('sort-by');
const resultsContainer = document.getElementById('results-container');
const resultCount = document.getElementById('result-count');
const loadingSpinner = document.getElementById('loading');
const articleTemplate = document.getElementById('article-template');

// 分页相关DOM
const paginationBox = document.getElementById('pagination');
const btnFirst   = document.getElementById('btn-first');
const btnPrev    = document.getElementById('btn-prev');
const btnNext    = document.getElementById('btn-next');
const btnLast    = document.getElementById('btn-last');
const pageInfo   = document.getElementById('page-info');
const pageJump   = document.getElementById('page-jump');
const btnJump    = document.getElementById('btn-jump');

// API URL
const API_BASE_URL = '/api';

// 分页全局变量
const ITEMS_PER_PAGE = 10;
let allArticles = [];   // 保存完整结果
let currentPage = 1;

// 初始化
document.addEventListener('DOMContentLoaded', function() {
    paginationBox.classList.add('hidden');
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

    // 分页按钮事件
    btnFirst.addEventListener('click', () => goToPage(1));
    btnPrev .addEventListener('click', () => goToPage(currentPage - 1));
    btnNext .addEventListener('click', () => goToPage(currentPage + 1));
    btnLast .addEventListener('click', () => goToPage(totalPages()));
    btnJump .addEventListener('click', () => {
        const p = parseInt(pageJump.value);
        if (!isNaN(p)) goToPage(p);
    });

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
    clearResults();
    paginationBox.classList.add('hidden');   // 先隐藏分页条

    // 清空结果容器
    clearResults();
    
    try {
        const params = new URLSearchParams({
            query: query,
            time_period: timePeriod.value,
            max_results: maxResults.value,
            sort: sortBy.value
        });
        const response = await fetch(`${API_BASE_URL}/articles?${params}`);
        const data = await response.json();

        if (data.status === 'success') {
            // 先按期刊分区过滤
            let filtered = data.data;
            if (journalRank.value !== 'all') {
                filtered = filterArticlesByRank(data.data, journalRank.value);
            }

            allArticles = filtered;        // 保存完整结果
            resultCount.textContent = `(${allArticles.length})`;
            currentPage = 1;               // 重置到第一页
            paginationBox.classList.add('hidden');
            if (allArticles.length === 0) {
                showEmptyState('未找到符合条件的文献');
            } else {
                renderPage();
                renderPagination();   // 里面会决定是否显示
            }
        } else {
            showEmptyState('搜索出错: ' + data.message);
        }
    } catch (error) {
        console.error('搜索出错:', error);
        showEmptyState('搜索请求失败，请稍后重试');
    } finally {
        showLoading(false);
    }
}

// 计算总页数
function totalPages() {
    return Math.ceil(allArticles.length / ITEMS_PER_PAGE) || 1;
}

// 跳转到指定页
function goToPage(page) {
    if (page < 1 || page > totalPages()) return;
    currentPage = page;
    renderPage();
    renderPagination();
    window.scrollTo({ top: 0, behavior: 'smooth' });
}

// 渲染当前页内容
function renderPage() {
    clearResults();
    if (allArticles.length === 0) {
        showEmptyState('未找到符合条件的文献');
        paginationBox.classList.add('hidden');
        return;
    }
    const start = (currentPage - 1) * ITEMS_PER_PAGE;
    const pageArticles = allArticles.slice(start, start + ITEMS_PER_PAGE);
    displayArticles(pageArticles);
}

// 渲染分页控件：只有结果条数 > ITEMS_PER_PAGE 才显示
function renderPagination() {
    // 结果不足一页时保持隐藏
    if (allArticles.length <= ITEMS_PER_PAGE) {
        paginationBox.classList.add('hidden');
        return;
    }
    paginationBox.classList.remove('hidden');

    pageInfo.textContent = `${currentPage} / ${totalPages()}`;
    btnFirst.disabled = currentPage === 1;
    btnPrev .disabled = currentPage === 1;
    btnNext .disabled = currentPage === totalPages();
    btnLast .disabled = currentPage === totalPages();
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
        const journalInfo = template.querySelector('.journal-info');
        
        // 创建第一行：期刊名称和分区
        const journalLine = document.createElement('div');
        journalLine.className = 'journal-line';
        
        const journalName = document.createElement('span');
        journalName.className = 'journal-name';
        journalName.textContent = article.journal.title || '期刊未知';
        journalLine.appendChild(journalName);
        
        // 添加分区信息（如果有）
        if (article.journal_metrics && article.journal_metrics.cas_division) {
            const casElement = document.createElement('span');
            casElement.className = 'cas-division';
            casElement.textContent = `中科院: ${article.journal_metrics.cas_division}`;
            journalLine.appendChild(casElement);
        }
        
        // 创建第二行：发表日期
        const publishDate = document.createElement('div');
        publishDate.className = 'publish-date';
        publishDate.textContent = article.journal.pub_date || '';
        
        // 清空原有内容并添加新的布局
        journalInfo.innerHTML = '';
        journalInfo.appendChild(journalLine);
        journalInfo.appendChild(publishDate);
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
       // 打印原始摘要内容到控制台，便于排查
        // console.log('原始摘要内容:', article.abstract);
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
    
    // 添加文章卡片点击事件，跳转到详情页面
    articleCard.style.cursor = 'pointer';
    articleCard.addEventListener('click', (e) => {
        // 确保点击的不是链接按钮
        if (!e.target.closest('.article-actions')) {
            window.open(`article-detail.html?pmid=${article.pmid}`, '_blank');
        }
    });
    
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