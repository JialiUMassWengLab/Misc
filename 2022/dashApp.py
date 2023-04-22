import re
import os
import sys
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from dash import Dash, Input, Output, dcc, html, no_update, ctx


data = pd.read_csv("data/STdeconvolve.K15.geneExpr.genes.TSNE.embedding.coor.csv",index_col=0)
celltypes = data.columns[2:].sort_values().unique()
data.loc[:,'geneName'] = data.index
genes = data['geneName'].sort_values().unique()
#image_path = "assets/dash.topicOnTissue.STdeconvolve.K15.geneExpr.None.png"
df_long = pd.melt(data.iloc[:,2:],id_vars='geneName',var_name='celltype')
data2 = pd.read_csv('data/snRNAseqData.csv',index_col=0)
data2 = data2[(data2['celltype']!='Unannotated') & (data2['celltype']!='Doublets')
              & (data2['celltype']!='Empty')]
data2.loc[:,'pct_detected'] = data2['pct_detected']*100

def get_blank_graph():
    return {
        "layout": {
            "xaxis": {
                "visible": False
            },
            "yaxis": {
                "visible": False
            },
            "annotations": [
                {
                    "text": "Please select a gene",
                    "xref": "paper",
                    "yref": "paper",
                    "showarrow": False,
                    "font": {"size": 28}
                }
            ]
        }
    }
    

external_stylesheets = [
    {
        "href": (
            "https://fonts.googleapis.com/css2?"
            "family=Lato:wght@400;700&display=swap"
        ),
        "rel": "stylesheet",
    },
]
app = Dash(__name__, external_stylesheets=external_stylesheets)
app.title = "Rat ovary Spatial Tx viewer"

app.layout = html.Div(
    children=[
        html.Div(
            children=[
                html.H1(
                    children="Rat ovary Spatial Tx gene viewer", className="header-title"
                ),
                html.P(
                    children=(
                        "Gene-centric visualization of rat ovary spatial tx data"
                    ),
                    className="header-description",
                ),
            ],
            className="header",
        ),
        html.Div(
            children=[
                html.Div(
                    children=[
                        html.Div(children="Inferred Cell Type", className="menu-title"),
                        dcc.Dropdown(
                            id="celltype-filter",
                            options=[
                                {"label": celltype, "value": celltype}
                                for celltype in celltypes
                            ],
                            #value="Adipocyte",
                            clearable=True,
                            searchable=False,
                            className="dropdown",
                        ),
                    ]
                ),
                html.Div(
                    children=[
                        html.Div(children="Gene", className="menu-title"),
                        dcc.Dropdown(
                            id="gene-filter",
                            options=[
                                {"label": gene, "value": gene}
                                for gene in genes
                            ],
                            clearable=True,
                            searchable=True,
                            className="dropdown",
                        ),
                    ]
                ),
            ],
            className="menu",
        ),
        html.Div(
            children=[
                html.Div(
                    id="tissue-image-celltype",
                    #children=html.Img(src=image_path),
                    className="card",
                ),
                html.Div(
                    children=dcc.Graph(
                        id="gene-scatterplot",
                        #config={"displayModeBar": False},
                        style={'height': '100vh'}
                    ),
                    style={'height': '100vh'},
                    className="card",
                ),
            ],
            className="wrapper",
            style={'width':'58vw','display':'inline-block'},
        ),
        html.Div(
            children=[
                html.Div(
                    children=dcc.Graph(
                        id="dist-bar-plot",
                        style={'height':'45vh'}
                    ),
                    className="card",
                    style={'height':'45vh'},
                ),
                html.Div(
                    children=dcc.Graph(
                        id="oralfsh-dot-plot",
                        style={'height':'50vh'}
                    ),
                    className="card",
                    style={'height':'50vh'},
                ),
                html.Div(
                    children=dcc.Graph(
                        id="optigon-dot-plot",
                        style={'height':'50vh'}
                    ),
                    className="card",
                    style={'height':'50vh'},
                ),
                html.Div(
                    id="tissue-image-gene",
                    className="card",
                ),
            ],
            className="wrapper",
            style={'width':'36vw','display':'inline-block',"verticalAlign": "top"},
        ),
    ]
)

@app.callback(
    Output("tissue-image-celltype", "children"),
    Output("gene-scatterplot", "figure"),
    Input("celltype-filter", "value"),
    Input("gene-filter","value"),
)
def update_celltype_info(celltype,gene):
    image_path = "assets/dash.topicOnTissue.STdeconvolve.K15.geneExpr.%s.png" % celltype
    tissue_image = html.Img(src=image_path,style={'height':'100%', 'width':'100%'})
    gene_scatterplot_figure = px.scatter(data,x='TSNE1',y='TSNE2',color=celltype,
                                         color_continuous_scale='reds',hover_name='geneName')
    if not gene==None:
        info = data.loc[gene,:]
        gene_scatterplot_figure.add_annotation(x=info['TSNE1'],y=info['TSNE2'],
                                               text=info['geneName'],arrowwidth=2,
                                               showarrow=True)

    return tissue_image, gene_scatterplot_figure


@app.callback(
    Output("dist-bar-plot", "figure"),
    Output("oralfsh-dot-plot", "figure"),
    Output("optigon-dot-plot", "figure"),
    Output("tissue-image-gene", "children"),
    Input("gene-scatterplot","hoverData"),
    Input("gene-filter","value"),
)
def update_gene_info(hoverData,gene_val):
    triggered_id = ctx.triggered_id
    gene = None
    
    if triggered_id == "gene-scatterplot":
        if not hoverData == None:
            gene = hoverData["points"][0]['hovertext']
    elif triggered_id == "gene-filter":
        gene = gene_val
        
    if gene == None:
        distribution_bar_plot = get_blank_graph()
        oralfsh_dot_plot = get_blank_graph()
        optigon_dot_plot = get_blank_graph()
        tissue_image2 = html.Br()
    else:
        df1 = df_long[df_long['geneName']==gene]
        distribution_bar_plot = px.bar(df1,x='celltype',y='value',title=gene)
        df2 = data2[data2['gene']==gene]
        oralfsh_df = df2[df2['study']=='oralFSH'].copy()
        optigon_df = df2[df2['study']=='Optigon'].copy()
        sizeref1 = max(oralfsh_df['pct_detected'])/max(df2['pct_detected'])*15
        sizeref2 = max(optigon_df['pct_detected'])/max(df2['pct_detected'])*15
        oralfsh_dot_plot = px.scatter(oralfsh_df,x='celltype',y='gene',size='pct_detected',
                                      color='avg_expr',size_max=sizeref1,
                                      color_continuous_scale='reds',title='oralFSH')
        optigon_dot_plot = px.scatter(optigon_df,x='celltype',y='gene',size='pct_detected',
                                      color='avg_expr',size_max=sizeref2,
                                      color_continuous_scale='reds',title='Optigon')
        
        image_path2 = "assets/dash.geneOnTissue.STdeconvolve.K15.geneExpr.%s.png" % gene
        tissue_image2 = html.Img(src=image_path2,style={'height':'100%','width':'100%'})
        
    return distribution_bar_plot, oralfsh_dot_plot, optigon_dot_plot, tissue_image2


if __name__ == "__main__":
    app.run_server("172.31.12.243",port=5015,debug=True)

