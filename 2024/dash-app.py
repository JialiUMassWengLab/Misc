import dash
import dash_design_kit as ddk
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output
from dash.dash_table.Format import Format
from dash import callback
import plotly.express as px
import pandas as pd
import re
import os
import boto3
from botocore.client import Config

theme = {
    "accent":"#3f4f75",
    "accent_positive":"#357e94",
    "accent_negative":"#f4564e",
    "background_content":"#f1f2f4",
    "background_page":"#fff",
    "body_text":"#20293D",
    "border":"white",
    "border_style":{
        "name":"underlined",
        "borderWidth":"0px 0px 1px 0px",
        "borderStyle":"solid",
        "borderRadius":0
    },
    "button_border":{
        "width":"1px",
        "color":"#3f4f75",
        "radius":"0px"
    },
    "button_capitalization":"uppercase",
    "button_text":"#3f4f75",
    "button_background_color":"#f1f2f4",
    "control_border":{
        "width":"0px 0px 1px 0px",
        "color":"white",
        "radius":"0px"
    },
    "control_background_color":"#f1f2f4",
    "control_text":"#20293D",
    "card_margin":"15px",
    "card_padding":"5px",
    "card_border":{
        "width":"0px 0px 0px 0px",
        "style":"solid",
        "color":"white",
        "radius":"0px"
    },
    "card_background_color":"#f1f2f4",
    "card_box_shadow":"0px 0px 0px rgba(0,0,0,0)",
    "card_outline":{
        "width":"0px",
        "style":"solid",
        "color":"white"
    },
    "card_header_margin":"0px",
    "card_header_padding":"10px",
    "card_header_border":{
        "width":"0px 0px 1px 0px",
        "style":"solid",
        "color":"white",
        "radius":"0px"
    },
    "card_header_background_color":"#f1f2f4",
    "card_header_box_shadow":"0px 0px 0px rgba(0,0,0,0)",
    "breakpoint_font":"1200px",
    "breakpoint_stack_blocks":"700px",
    "colorway":[
        "#3f4f75",
        "#80cfbe",
        "#60bab6",
        "#f4564e",
        "#ffeeb2",
        "#20293d",
        "#faddd2",
        "#ffdd68",
        "#357e94",
        "#a1acc3"
    ],
    "colorscale":[
        "#ffffff",
        "#f0f0f0",
        "#d9d9d9",
        "#bdbdbd",
        "#969696",
        "#737373",
        "#525252",
        "#252525",
        "#000000"
    ],
    "dbc_primary":"#3f4f75",
    "dbc_secondary":"#a6a6a6",
    "dbc_info":"#00769F",
    "dbc_gray":"#adb5bd",
    "dbc_success":"#00CCA4",
    "dbc_warning":"#FADD6A",
    "dbc_danger":"#F76065",
    "font_family":"Open Sans",
    "font_family_header":"Domine",
    "font_family_headings":"Domine",
    "font_size":"17px",
    "font_size_smaller_screen":"15px",
    "font_size_header":"24px",
    "title_capitalization":"capitalize",
    "header_content_alignment":"spread",
    "header_margin":"0px 0px 15px 0px",
    "header_padding":"0px",
    "header_border":{
        "width":"0px 0px 0px 0px",
        "style":"solid",
        "color":"white",
        "radius":"0px"
    },
    "header_background_color":"#f1f2f4",
    "header_box_shadow":"none",
    "header_text":"#20293D",
    "heading_text":"#20293D",
    "text":"#20293D",
    "report_background_content":"#FAFBFC",
    "report_background_page":"white",
    "report_text":"black",
    "report_font_family":"Computer Modern",
    "report_font_size":"12px"
}

gwas_summary_row_colors = [
                {
                    'if': {
                        'filter_query': '{{GWAS evidence}} = {}'.format('Strong'),
                    },
                    'backgroundColor': '#229954',
                    'color': 'white'
                },
                {
                    'if': {
                        'filter_query': '{{GWAS evidence}} = {}'.format('Strong/Medium'),
                    },
                    'backgroundColor': '#ABEBC6',
                },
                {
                    'if': {
                        'filter_query': '{{GWAS evidence}} = {}'.format('Medium'),
                    },
                    'backgroundColor': '#FAD7A0',
                },
                {
                    'if': {
                        'filter_query': '{{GWAS evidence}} = {}'.format('Medium/Low'),
                    },
                    'backgroundColor': '#DC7633',
                    'color': 'white'
                },
                {
                    'if': {
                        'filter_query': '{{GWAS evidence}} = {}'.format('Low'),
                    },
                    'backgroundColor': '#A04000',
                    'color': 'white'
                },
            ]


def appendTab(df,tableName,tabList):
    tab = dcc.Tab(label=tableName,children=[
            ddk.Card(width = 100, children=[
                ddk.CardHeader(title=tableName),
                ddk.DataTable(
                    id='%s table' % tableName,
                    #columns=[{"name": i, "id": i} if not i=='External Link' else {'id': i, 'name': i, 'presentation': 'markdown'} for i in df.columns[1:]],
                    columns=[{"name": i, "id": i} if not i=='External Link' and df.dtypes[i]==object \
                        else {"name": i, "id": i, "type":'numeric', "format":Format()} if not i=='External Link' \
                            else {'id': i, 'name': i, 'presentation': 'markdown'} for i in df.columns],
                    data=df.to_dict('records'),
                    style_header={
                        'whiteSpace': 'normal',
                        'height': 'auto',
                        'fontWeight': 'bold'
                    },
                    style_table={'overflowX': 'auto'},
                    style_cell={'padding': '5px'},
                    page_size=20,
                    filter_action="native",
                    sort_action="native",
                    editable=False
                )
            ])
        ], id='%s-tab' % tableName)

    return tabList.append(tab)


# load your AWS S3 credentials from environment variables saved to your dash app
access_key = os.environ.get("S3_ACCESS_KEY")
secret_key = os.environ.get("S3_SECRET_KEY")
bucket_name = os.environ.get("S3_BUCKET_NAME")
path_prefix = 'home/chris_deboever/gsp_ie/20231203/'
# connect to S3 using boto
client = boto3.client(
    "s3",
    aws_access_key_id=access_key,
    aws_secret_access_key=secret_key,
    config=Config(signature_version="s3v4"),
)

# get a list of all genes
#resp = client.list_objects_v2(Bucket=bucket_name,Prefix=path_refix)
#geneList = [os.path.basename(x['Key']).replace(' GSP.xlsx','') for x in resp.get('Contents') if '.xlsx' in x['Key']]
paginator = client.get_paginator('list_objects_v2')
pages = paginator.paginate(Bucket=bucket_name, Prefix=path_prefix)
geneList = []
for page in pages:
    geneList += [os.path.basename(x['Key']).replace(' GSP.xlsx','') for x in page['Contents'] if '.xlsx' in x['Key']]
#print(geneList)

#get gene info
array = []
for gene in geneList:
    gene_path = os.path.join(path_prefix,gene+' GSP','GOI Gene Info.tsv.gz')
    gene_df = pd.read_csv(client.get_object(Bucket=bucket_name, Key=gene_path)['Body'],sep='\t',compression='gzip')
    array.append(gene_df)
GOI_df = pd.concat(array,join='inner',axis=0)

#get sheet description
desc_path = os.path.join(path_prefix,geneList[0]+' GSP','Sheet Descriptions.tsv.gz')
desc = pd.read_csv(client.get_object(Bucket=bucket_name, Key=desc_path)['Body'],sep='\t',index_col=0,compression='gzip')
tooltips = [dbc.Tooltip(row[1],target='%s-tab' % row[0]) for row in desc.iterrows()]

#read GWAS summary dataframe
gwas_path = os.path.join(path_prefix,'aggregated_GWAS_summary.csv.gz')
summary_df = pd.read_csv(client.get_object(Bucket=bucket_name, Key=gwas_path)['Body'],compression='gzip')
#generate dataframes for figures
counts = summary_df.groupby(['GOI','GWAS evidence']).size().reset_index(name='counts')
counts_array = []
for pair,group in summary_df.groupby(['GOI','Trait category']):
    counts_array.append(group.sort_values('Max GOI Evidence Strength rank',ascending=False).iloc[0,:])
counts2 = pd.concat(counts_array,join='inner',axis=1).transpose()
#print(counts2)
#make figures
fig_bar = px.bar(counts,x='GOI',y='counts',color='GWAS evidence',title='GWAS supported traits',\
    labels={'counts':'Number of traits','GOI':'Gene of interest'},\
    color_discrete_map={
        'Strong': '#229954',
        'Strong/Medium': '#ABEBC6',
        'Medium': '#FAD7A0',
        'Medium/Low': '#DC7633',
        'Low': '#A04000'},\
    category_orders={'GWAS evidence':['Strong','Strong/Medium','Medium','Medium/Low','Low']})
fig_bar.update_layout(xaxis_tickangle=-45, xaxis={'categoryorder':'total descending'})
fig_dot = px.scatter(counts2,x='GOI',y='Trait category',color='GWAS evidence',\
    title='Trait categories supported', labels={'GOI':'Gene of interest'},
    color_discrete_map={
        'Strong': '#229954',
        'Strong/Medium': '#ABEBC6',
        'Medium': '#FAD7A0',
        'Medium/Low': '#DC7633',
        'Low': '#A04000'},\
    category_orders={'GWAS evidence':['Strong','Strong/Medium','Medium','Medium/Low','Low']})
fig_dot.update_layout(xaxis_tickangle=-45)

app = dash.Dash(__name__)
server = app.server  # expose server variable for Procfile

app.layout = ddk.App(show_editor=False, children=[
    ddk.Header([ddk.Title('Genetic Support Profiler')]),
    ddk.Card(width=95,children=[
        ddk.Graph(figure=fig_bar),
        ddk.Graph(figure=fig_dot),
    ]),
    ddk.ControlCard(width=95,children=[
        ddk.ControlItem(
            dcc.Dropdown(
                id='selected-gene-of-interest',
                options=[
                    {'label':i, 'value':i} for i in geneList
                ],
                multi=False,
                searchable=True
            ),
            label='Gene of Interest'
        )
    ]),
    ddk.Card(id='gene-info-table', width=95),
    ddk.Card(width=95, children=[
        ddk.CardHeader(title='GWAS Support Summary'),
        ddk.ControlCard(width=40,children=[
            ddk.ControlItem(
                dcc.Dropdown(
                    id='selected-trait-category',
                    options=[
                        {'label':i, 'value':i} for i in summary_df['Trait category'].unique()
                    ],
                    multi=False,
                    searchable=True,
                    placeholder='Filter by trait category'
                ),
            ),
        ]),
        ddk.DataTable(
            id='evidence support table',
            tooltip_header = {'GWAS evidence': 'strong: variant in credible set is a protein-coding variant in gene of interest\n \
                strong/medium: lead variant or poxy is protein-coding in gene of interest (fine-mapping not available)\n \
                medium: colocalization/credible set overlap with relevant tissue eQTL'},
            style_table={'overflowX':'auto'},
            style_cell={'padding': '5px'},
            style_header={
                'whiteSpace': 'normal',
                'height': 'auto',
                'fontWeight': 'bold'
            },
            page_size=30,
            filter_action="native",
            sort_action="native",
            style_data_conditional=gwas_summary_row_colors,
            editable=False
        ),
        html.Br(),
        html.Details([
            html.Summary('GWAS related data tables',style={'font-weight': 'bold'}),
            html.Br(),
            dcc.Tabs(id='gwas-related-tabs'),
        ])
    ]),
    ddk.Card(width=95, children=[
        ddk.CardHeader(title='OMIM Support Summary'),
        ddk.DataTable(
            id='omim summary table',
            style_table={'overflowX':'auto'},
            style_cell={'padding': '5px', 'textAlign': 'left'},
            style_header={
                'whiteSpace': 'normal',
                'height': 'auto',
                'fontWeight': 'bold'
            },
            page_size=30,
            filter_action="native",
            sort_action="native",
            editable=False
        ),
        html.Br(),
        html.Details([
            html.Summary('Mendelian diseases data tables',style={'font-weight': 'bold'}),
            html.Br(),
            dcc.Tabs(id='mendelian-disease-tabs'),
        ])
    ]),
    ddk.Card(width=95, children=[
        ddk.CardHeader(title='Gene burden Support Summary'),
        ddk.DataTable(
            id='rvas summary table',
            style_table={'overflowX':'auto'},
            style_cell={'padding': '5px'},
            style_header={
                'whiteSpace': 'normal',
                'height': 'auto',
                'fontWeight': 'bold'
            },
            page_size=30,
            filter_action="native",
            sort_action="native",
            editable=False
        ),
        html.Br(),
        html.Details([
            html.Summary('Gene burden and exome variant associations data tables',style={'font-weight': 'bold'}),
            html.Br(),
            dcc.Tabs(id='gene-burden-tabs'),
        ])
    ]),
    ddk.Card(width=95, children=[
        html.Details([
            html.Summary('Genetic instrument data tables',style={'font-weight': 'bold'}),
            html.Br(),
            dcc.Tabs(id='genetic-instrument-tabs'),
        ])
    ]),
    ddk.Card(width=95, children=[
        html.Details([
            html.Summary('Mouse phenotypes data tables',style={'font-weight': 'bold'}),
            html.Br(),
            dcc.Tabs(id='mouse-phenotype-tabs'),
        ])
    ]),
    html.Div(tooltips)
], theme=theme)


@callback(
    Output('gene-info-table','children'),
    Output('evidence support table','data'),
    Output('evidence support table','columns'),
    Output('omim summary table','data'),
    Output('omim summary table','columns'),
    Output('rvas summary table','data'),
    Output('rvas summary table','columns'),
    Output('gwas-related-tabs','children'),
    Output('mendelian-disease-tabs','children'),
    Output('gene-burden-tabs','children'),
    Output('genetic-instrument-tabs','children'),
    Output('mouse-phenotype-tabs','children'),
    Input('selected-gene-of-interest','value'),
    Input('selected-trait-category','value')
)
def render_tabs(GOI,trait):
    #gene info table
    tmp_df = GOI_df[GOI_df['Gene Symbol']==GOI].copy()
    geneInfoChildren = [
        ddk.CardHeader(title='Gene information'),
        ddk.DataTable(
            id='gene table',
            columns=[{"name": i, "id":i} for i in tmp_df.columns[1:12]],
            data=tmp_df.iloc[:,1:12].to_dict('records'),
            style_table={'overflowX':'auto'},
            style_cell={'padding': '5px'},
            style_header={
                'whiteSpace': 'normal',
                'height': 'auto',
                'fontWeight': 'bold'
            },
            editable=False
        )
    ]

    #GWAS support summary table
    tmp_df2 = pd.DataFrame(columns=['Trait category','Trait Reported','GWAS evidence'])
    if not GOI==None:
        table_path = os.path.join(path_prefix,GOI+' GSP','GWAS Topline Summary.tsv.gz')
        tmp_df2 = pd.read_csv(client.get_object(Bucket=bucket_name, Key=table_path)['Body'],sep='\t',compression='gzip')
        if not trait==None:
            tmp_df2 = tmp_df2[tmp_df2['Trait category']==trait]
    
    #OMIM summary table
    omim_summary = pd.DataFrame(columns=["OMIM ID","Disease Name","Category"])
    if not GOI==None:
        table_path = os.path.join(path_prefix,GOI+' GSP','Mendelian Topline Summary.tsv.gz')
        omim_summary = pd.read_csv(client.get_object(Bucket=bucket_name, Key=table_path)['Body'],sep='\t',compression='gzip')

    #RVAS summary table
    rvas_summary = pd.DataFrame(columns=['Trait category','Trait reported','Min logP'])
    if not GOI==None:
        table_path = os.path.join(path_prefix,GOI+' GSP','RVAS Topline Summary.tsv.gz')
        rvas_summary = pd.read_csv(client.get_object(Bucket=bucket_name, Key=table_path)['Body'],sep='\t',compression='gzip')

    #detailed tables in tabs
    tabList = []
    for tableName in ['GWAS Summary','OTG GWAS Assocs','OTG Lead Var Annot','GWAS Credible Sets','LD Neighbors','Missense Annotations','VEP transcript']:
        if not GOI==None:
            table_path = os.path.join(path_prefix,GOI+' GSP',tableName + '.tsv.gz')
            try:
                df = pd.read_csv(client.get_object(Bucket=bucket_name, Key=table_path)['Body'],sep='\t',compression='gzip')
            except botocore.exceptions.ClientError as e:
                print('Error %s when loading table %s' % (e.response['Error']['Code'],e.response['Error']['Key']))
        else:
            table_path = os.path.join(path_prefix,geneList[0]+' GSP',tableName + '.tsv.gz')
            df0 = pd.read_csv(client.get_object(Bucket=bucket_name, Key=table_path)['Body'],sep='\t',compression='gzip')
            df = pd.DataFrame(columns=df0.columns)
        appendTab(df,tableName,tabList)

    tabList2 = []
    for tableName in ['OMIM','OMIM HPO']:
        if not GOI==None:
            table_path = os.path.join(path_prefix,GOI+' GSP',tableName + '.tsv.gz')
            try:
                df = pd.read_csv(client.get_object(Bucket=bucket_name, Key=table_path)['Body'],sep='\t',compression='gzip')
            except botocore.exceptions.ClientError as e:
                print('Error %s when loading table %s' % (e.response['Error']['Code'],e.response['Error']['Key']))
        else:
            table_path = os.path.join(path_prefix,geneList[0]+' GSP',tableName + '.tsv.gz')
            df0 = pd.read_csv(client.get_object(Bucket=bucket_name, Key=table_path)['Body'],sep='\t',compression='gzip')
            df = pd.DataFrame(columns=df0.columns)
        appendTab(df,tableName,tabList2)

    tabList3 = []
    for tableName in ['UKB QT Burden','UKB Binary Burden','UKB QT Variant','UKB Binary Variant','OT Gene Burden']:
        if not GOI==None:
            table_path = os.path.join(path_prefix,GOI+' GSP',tableName + '.tsv.gz')
            try:
                df = pd.read_csv(client.get_object(Bucket=bucket_name, Key=table_path)['Body'],sep='\t',compression='gzip')
            except botocore.exceptions.ClientError as e:
                print('Error %s when loading table %s' % (e.response['Error']['Code'],e.response['Error']['Key']))
        else:
            table_path = os.path.join(path_prefix,geneList[0]+' GSP',tableName + '.tsv.gz')
            df0 = pd.read_csv(client.get_object(Bucket=bucket_name, Key=table_path)['Body'],sep='\t',compression='gzip')
            df = pd.DataFrame(columns=df0.columns)
        appendTab(df,tableName,tabList3)

    tabList4 = []
    for tableName in ['eQTL Variants','pQTL Variants','Predicted Functional Variants']:
        if not GOI==None:
            table_path = os.path.join(path_prefix,GOI+' GSP',tableName + '.tsv.gz')
            try:
                df = pd.read_csv(client.get_object(Bucket=bucket_name, Key=table_path)['Body'],sep='\t',compression='gzip')
            except botocore.exceptions.ClientError as e:
                print('Error %s when loading table %s' % (e.response['Error']['Code'],e.response['Error']['Key']))
        else:
            table_path = os.path.join(path_prefix,geneList[0]+' GSP',tableName + '.tsv.gz')
            df0 = pd.read_csv(client.get_object(Bucket=bucket_name, Key=table_path)['Body'],sep='\t',compression='gzip')
            df = pd.DataFrame(columns=df0.columns)
        appendTab(df,tableName,tabList4)

    tabList5 = []
    for tableName in ['Mouse Model Phenos']:
        if not GOI==None:
            table_path = os.path.join(path_prefix,GOI+' GSP',tableName + '.tsv.gz')
            try:
                df = pd.read_csv(client.get_object(Bucket=bucket_name, Key=table_path)['Body'],sep='\t',compression='gzip')
            except botocore.exceptions.ClientError as e:
                print('Error %s when loading table %s' % (e.response['Error']['Code'],e.response['Error']['Key']))
        else:
            table_path = os.path.join(path_prefix,geneList[0]+' GSP',tableName + '.tsv.gz')
            df0 = pd.read_csv(client.get_object(Bucket=bucket_name, Key=table_path)['Body'],sep='\t',compression='gzip')
            df = pd.DataFrame(columns=df0.columns)
        appendTab(df,tableName,tabList5)

    return geneInfoChildren, tmp_df2.to_dict('records'), [{"name": i, "id":i} for i in tmp_df2.columns], \
        omim_summary.to_dict('records'), [{"name": i, "id":i} for i in omim_summary.columns], \
            rvas_summary.to_dict('records'), [{"name": i, "id":i} for i in rvas_summary.columns], \
                tabList, tabList2, tabList3, tabList4, tabList5


if __name__ == '__main__':
    app.run_server(debug=True)
