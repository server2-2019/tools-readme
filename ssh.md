#SSH command
=====

## Table of Contents
<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
- [ssh-keyscan](#ssh-keyscan)
- [taxonkit](#taxonkit)
- [list](#list)
<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## ssh-keyscan
Source:
- [ssh-keyscan 命令](https://blog.csdn.net/airfish2000/article/details/55654172)<br>
```
ssh-keyscan [选项] [主机|地址列表名称列表]
```
| 选项  | 含义 |
| ------------- | ------------- |
| -H  | 在输出中哈希所有主机名和地址  |
| -p <端口号>  | 指定连接远程主机的端口号  |
| -T<超时时间>  | 设置连接尝试的超时时间  |
| -v  | 显示详细信息  |
| -t <类型>  | 指定从扫描的主机获取密钥的类型。类型可以是rsa或dsa，默认值是rsa  |

Examples
### ssh-keyscan -t rsa 192.30.253.113
```
# 192.30.253.113:22 SSH-2.0-babeld-959a4830
192.30.253.113 ssh-rsa AAAAB3NzaC1yc2EAAAABIwAAAQEAq2A7hRGmdnm9tUDbO9IDSwBK6TbQa+PXYPCPy6rbTrTtw7PHkccKrpp0yVhp5HdEIcKr6pLlVDBfOLX9QUsyCOV0wzfjIJNlGEYsdlLJizHhbn2mUjvSAHQqZETYP81eFzLQNnPHt4EVVUh7VfDESU84KezmD5QlWpXLmvU31/yMf+Se8xhHTvKSCZIFImWwoG6mbUoWf9nzpIoaSjB+weqqUUmpaaasXVal72J+UX2B+2RPW3RcT0eOzQgqlJL3RKrTJvdsjE3JEAvGq3lGHSZXy28G3skua2SmVi/w4yCE6gbODqnTWlg7+wC604ydGXA8VJiS5ap43JXiUFFAaQ==
```
***
## taxonkit
TaxonKit - A Cross-platform and Efficient NCBI Taxonomy Toolkit<br>
Version: 0.5.0<br>
```
```
Examples
### 1. Default usage
```
$ taxonkit list --ids 9605,239934
  ```
