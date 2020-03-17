#SSH command
=====

## Table of Contents
<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
- [ssh-keyscan](#ssh-keyscan)
- [ssh-keygen](#ssh-keygen)
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
## ssh-keygen
Source:
- [SSH key的介绍与在Git中的使用](https://www.jianshu.com/p/1246cfdbe460)<br>
ssh-keygen命令用于为ssh生成、管理和转换认证密钥，它支持RSA和DSA两种认证密钥。<br>
该命令的选项：
```
-b：指定密钥长度；
-e：读取openssh的私钥或者公钥文件；
-C：添加注释；
-f：指定用来保存密钥的文件名；
-i：读取未加密的ssh-v2兼容的私钥/公钥文件，然后在标准输出设备上显示openssh兼容的私钥/公钥；
-l：显示公钥文件的指纹数据；
-N：提供一个新密语；
-P：提供（旧）密语；
-q：静默模式；
-t：指定要创建的密钥类型。

生成密钥对时，有一个选项要求你设置密码（passphrase），该密码是用来保护你的私钥的密码。
如果设置了则在使用私钥时会要求你输入这个密码；一般不设置，记不住【之后还可更改此密码，使用ssh-keygen -p】。

生成后最好将私钥进行备份。另还有-C选项，用于为指定注释，通常使用自己的邮件名作为注释。
-b bits选项 Specifies the number of bits in the key to create. For RSA keys, the minimum size 
is 1024 bits and the default is 2048 bits. Generally, 2048 bits is considered sufficient.
DSA keys must be exactly 1024 bits

作者：faner
链接：https://www.jianshu.com/p/1246cfdbe460
来源：简书
著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。
```
Examples
### 1. 示例：为了安全考虑使用RSA加密方式并指定密钥长度 -b 2048（1024的密钥长度能够被破解，建议指定为2048或4096）
```
$ ssh-keygen -t rsa -C "your_email@example.com" -b 2048
Generating public/private rsa key pair.
Enter file in which to save the key
(/Users/your_user_directory/.ssh/id_rsa): 按回车键 （如果需要生成多对key，则输入/home/users/.ssh/filename）
Enter passphrase (empty for no passphrase): 输入密码(一般不输入密码，直接回车)
Enter same passphrase again: 再次输入密码
...

# 查看公钥文件中的内容
$ cat ~/.ssh/id_rsa.pub
ssh-rsa "公钥内容" your_email@example.com

# 注意在其他地方导入公钥时一定要将公钥文件中的*全部内容*都导入，包括末尾你的邮箱。
```
### 2. 实际操作的一次示例：
```
$ ssh-keygen -t rsa -C "Fan@outlook.com" -b 2048
Generating public/private rsa key pair.
Enter file in which to save the key (/home/fan/.ssh/id_rsa): /home/fan/.ssh/FDGitHub_rsa
Enter passphrase (empty for no passphrase): 
Enter same passphrase again: 
Your identification has been saved in /home/fan/.ssh/FDGitHub_rsa.
Your public key has been saved in /home/fan/.ssh/FDGitHub_rsa.pub.
The key fingerprint is:
SHA256:GcK7ORvFzH6fzA7qPmnzBr1DOWho5cCVgIpLkh6VGb8 Fan@outlook.com
The key's randomart image is:
+---[RSA 2048]----+
|   .+... .       |
|   +o.  o        |
| o.. oo..        |
|+o.   +*.o       |
|+..  E.=So .     |
|..    o== =      |
|     .=..+oo     |
|       +=o+= .   |
|      .++=.o*    |
+----[SHA256]-----+
```
### 3. GitHub使用非默认密钥对:
```
Host github.com
RSAAuthentication yes
# 也可以使用公钥
IdentityFile ~/.ssh/FDGitHub_rsa.pub
```
配置完成后可以使用如下命令测试连接：
```
# 测试时替换掉 example.com
ssh -T git@example.com
# 例如 gitlab
ssh -T git@gitlab.com
# 例如 github
ssh -T git@github.com
# 例如 coding
ssh -T git@git.coding.net
# 例如 码云
ssh -T git@gitee.com
# bitbucket
ssh -T git@bitbucket.org

# 也可以使用下面的命令来调试连接
ssh -Tv git@example.com
```
