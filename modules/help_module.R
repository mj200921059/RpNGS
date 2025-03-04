# Help Tab 模块

# UI 部分
mod_help_ui <- function(id) {
  ns <- NS(id)  # 用于命名空间，防止 ID 冲突
  tabPanel(h4('Help'), value = 'help',
           shiny::tags$style(
             '.navbar-nav {width: 90%;} .navbar-nav :first-child{float:right}'
           ),
           column(1),
           column(10,
                  div(style = "color:#337ab7; text-align:center;", titlePanel(h2("RpNGS Operation Details"))),
                  br(),
                  # 添加滚动区域
                  shiny::tags$style(HTML("
                        .markdown-container {
                        width: 100%;
                        height: 900px;
                        overflow: auto;
                        border: 1px solid #ccc;
                        }
                  ")),
                  shiny::tags$div(class = "markdown-container",
                                  includeMarkdown(paste0(getwd(),"/www/html/help.Rmd"))  # Markdown 文件路径
                  )
           ),
           column(1)
  )
}

# 服务器逻辑 (可以为空，如果 Help 页不需要交互)
mod_help_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    # 暂时没有需要处理的交互逻辑
  })
}
