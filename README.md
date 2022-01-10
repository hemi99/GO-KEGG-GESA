##### 数据来源
1.数据来源：这是一个对于当前热门的癌症有关的基因，例如TP53等。
##### 富集分析的三种类型
2.富集分析：包括了GO富集分析、KEGG富集分析以及GESA富集分析
##### GO与KEGG富集分析结果
3.富集分析结果：富集分析每个部分得到的结果存放于xls文件中
##### 富集分析意义
4.意义：富集到的terms有助于后续的研究分析，可以通过表型来验证基因富集到的terms的准确性
###### 富集结果输出文件说明
>>all_enrich.xls文件：存放的是取函数enrichGO其参数为ont="all"时的富集结果

>>BP_enrich.xls文件：存放的是取函数enrichGO其参数为ont="BP"时的富集结果

>>MF_enrich.xls文件：存放的是取函数enrichGO其参数为ont=""时的富集结果

>>CC_enrich.xls文件：存放的是取函数enrichGO其参数为ont="MF"时的富集结果

>>KEGG_enrich.xls文件：存放的是函数enrichKEGG做KEGG时的富集结果
