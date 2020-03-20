#Python functions

=====

## Table of Contents
<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
- [判断文件是否存在](#判断文件是否存在)
- [去字符串收尾空格](#去字符串收尾空格)
- [判断两个字符串是否相同](#判断两个字符串是否相同)
- [python值相同变量不同内存值是否相同](#python值相同变量不同内存值是否相同)
- [字符串拼接](#字符串拼接)
<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Python
Source:
- [Python3.8 标准库](https://docs.python.org/zh-cn/3.8/library/)<br>
- [Python PEP](https://www.python.org/dev/peps/)<br>
- [Python中下划线的5种含义](https://zhuanlan.zhihu.com/p/36173202)
- 

```
xxx
```
| 选项  | 含义 |
| ------------- | ------------- |
| xxx  | xxx  |

***
## ssh-keygen
Source:
- [判断文件是否存在的三种做法](https://www.cnblogs.com/jhao/p/7243043.html)<br>

```

```
Examples
### 1. XXX
```
XXX
```
***
## 去字符串收尾空格
Source:
- [去字符串收尾空格](https://www.jianshu.com/p/bd953fde69e6)<br>

```
字符串“  illusion  ”首尾都有空格
去左边空格，用lstrip()
去右边空格，用rstrip()
去首尾空格，用strip()
```
## 判断两个字符串是否相同
Source:
- [判断两个字符串是否相同](https://blog.csdn.net/weixin_34146805/article/details/85828509)<br>
```
The operators ``is`` and ``is not`` test for object identity: ``x is
y`` is true if and only if *x* and *y* are the same object.  ``x is
not y`` yields the inverse truth value.
```
### 1. is 主要是判断 2 个变量是否引用的是同一个对象，如果是的话，则返回 true，否则返回 false。
```
 is 相等代表两个对象的 id 相同（从底层来看的话，可以看作引用同一块内存区域）
```
Examples
### 2. == 用来判断两个对象的值是否相等
Examples
```
#-*-conding:utf-8-*-
i='xinwen';
m=input();
if i==m:
    print('yes');
else:
    print('no');
 
input();
```
***
## python值相同变量不同内存值是否相同
Source:
- [python 值相同变量名不同，内存地址相同吗？](https://bbs.csdn.net/topics/392267652)<br>
小整数对象池：python在执行的时候，为了节约空间，帮我们创建好了小整数对象池，[-5~256]，都是固定的地址，不管你用不用，都会存在。
```
比如，a=5,b=5,id(a)和id(b)的地址都是小整数池中固定已经存在的地址，所以相等
但如果，a=1000,b=1000,这时会给a一个地址，也会给b一个地址，但他们都不相等。
```
Examples
```
>>> a = 256
>>> b = 256
>>> id(a)
9987148
>>> id(b)
9987148
>>> a = 257
>>> b = 257
>>> id(a)
11662816
>>> id(b)
11662828
```
***
## 字符串拼接
Source:
- [python字符串拼接方法](https://blog.csdn.net/weixin_39274753/article/details/81660887)<br>
### 1. 方法1 “+” 拼接
Examples
```
name = 'Jack'
age = 18
merge = name + "\t" + age
```
### 2. 方法2 “%” 拼接
Examples
```
name = 'Jack'
age = 18
merge = "%s\t%s" % (name,age)
```
### 3. 方法3 “format” 拼接
Examples
```
name = 'Jack'
age = 18
merge = "{}\t{}".format(name,age)
# or
merge = "{0}\t{1}".format(name,age)
```
### 4. 方法4 “join” 拼接
Examples
```
name = 'Jack'
age = 18
merge = "\t".join([name,str(age)])
```
***
