#Python functions

=====

## Table of Contents
<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
- [object class](#object-class)
- [argparse CustomAction](#argparse-customaction)
- [判断文件是否存在](#判断文件是否存在)
- [去字符串收尾空格](#去字符串收尾空格)
- [判断两个字符串是否相同](#判断两个字符串是否相同)
- [python值相同变量不同内存值是否相同](#python值相同变量不同内存值是否相同)
- [字符串拼接](#字符串拼接)
- [列表操作](#列表操作)
- [字典操作](#字典操作)
- [字符串操作](#字符串操作)
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
## object class
- [object class](https://blog.csdn.net/qq_34964399/article/details/79119896)
```
class Cat:
  def __init__(self):
    print("-----haha-----")

  def eat(self):
    print("cat is eating a fish ...")

  def drink(self):
    print("cat is drinking kele ...")

  def introduce(self):
    print("%s is %d years old" % (self.name, self.age))

# create a object
tom = Cat()
tom.eat()
tom.drink()
tom.name = 'Tom'
tom.age = 40
tom.introduce()

lanmao = Cat()
lanmao.name = "blue cat"
lanmao.age = 10
lanmao.introduce()

```
## argparse CustomAction
Source:
- [argparse CustomAction](https://pymotw.com/2/argparse/)
```
import argparse

class CustomAction(argparse.Action):
    def __init__(self,
                 option_strings,
                 dest,
                 nargs=None,
                 const=None,
                 default=None,
                 type=None,
                 choices=None,
                 required=False,
                 help=None,
                 metavar=None):
        argparse.Action.__init__(self,
                                 option_strings=option_strings,
                                 dest=dest,
                                 nargs=nargs,
                                 const=const,
                                 default=default,
                                 type=type,
                                 choices=choices,
                                 required=required,
                                 help=help,
                                 metavar=metavar,
                                 )
        print
        print 'Initializing CustomAction'
        for name,value in sorted(locals().items()):
            if name == 'self' or value is None:
                continue
            print '  %s = %r' % (name, value)
        return

    def __call__(self, parser, namespace, values, option_string=None):
        print
        print 'Processing CustomAction for "%s"' % self.dest
        print '  parser = %s' % id(parser)
        print '  values = %r' % values
        print '  option_string = %r' % option_string
        
        # Do some arbitrary processing of the input values
        if isinstance(values, list):
            values = [ v.upper() for v in values ]
        else:
            values = values.upper()
        # Save the results in the namespace using the destination
        # variable given to our constructor.
        setattr(namespace, self.dest, values)

parser = argparse.ArgumentParser()

parser.add_argument('-a', action=CustomAction)
parser.add_argument('-m', nargs='*', action=CustomAction)
parser.add_argument('positional', action=CustomAction)

results = parser.parse_args(['-a', 'value', '-m' 'multi-value', 'positional-value'])
print
print results
```

## 去字符串收尾空格
Source:
- [去字符串收尾空格](https://www.jianshu.com/p/bd953fde69e6)<br>

```
字符串“  illusion  ”首尾都有空格
去左边空格，用lstrip()
去右边空格，用rstrip()
去首尾空格，用strip()
```
***
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
## 列表操作
Source:
- [Python 列表(List)操作方法详解](https://blog.csdn.net/zhu_liangwei/article/details/7931701)<br>
### 1. 创建列表
只要把逗号分隔的不同的数据项使用方括号括起来即可。如下所示：
```
list1 = ['physics', 'chemistry', 1997, 2000];
list2 = [1, 2, 3, 4, 5 ];
sample_list = ['a',1,('a','b')]
```
### 2. 访问列表
使用下标索引来访问列表，截取列表，删除列表的第一个值，列表中插入一个值，遍历列表，如下所示：
```
sample_list = ['a','b',0,1,3]
#访问列表
value_start = sample_list[0]
end_value = sample_list[-1]
# 截取列表
sample_list[1:5]
# 删除列表的第一个值
del sample_list[0]
# 列表中插入一个值
sample_list[0:0] = ['sample value']
# 得到列表的长度
list_length = len(sample_list)
# 列表遍历
for element in sample_list:
  print(element)
```
说明
```
以0开始，有负下标的使用
0第一个元素，-1最后一个元素，
-len第一个元 素，len-1最后一个元素
取list的元素数量
len(list) #list的长度。实际该方法是调用了此对象的__len__(self)方法。
list的方法
L.append(var) #追加元素
L.insert(index,var)
L.pop(var) #返回最后一个元素，并从list中删除之
L.remove(var) #删除第一次出现的该元素
L.count(var) #该元素在列表中出现的个数
L.index(var) #该元素的位置,无则抛异常
L.extend(list) #追加list，即合并list到L上
L.sort() #排序
L.reverse() #倒序
list 操作符:,+,*，关键字del
a[1:] #片段操作符，用于子list的提取
[1,2]+[3,4] #为[1,2,3,4]。同extend()
[2]*4 #为[2,2,2,2]
del L[1] #删除指定下标的元素
del L[1:3] #删除指定下标范围的元素
list的复制
L1 = L #L1为L的别名，用C来说就是指针地址相同，对L1操作即对L操作。函数参数就是这样传递的
L1 = L[:] #L1为L的克隆，即另一个拷贝。

```
***
## 字典操作
Source:
- [字典操作](https://blog.csdn.net/zhu_liangwei/article/details/7931701)<br>
### 1. 常规操作

```
dict = {‘ob1′:’computer’, ‘ob2′:’mouse’, ‘ob3′:’printer’}
每一个元素是pair，包含key、value两部分。key是Integer或string类型，value 是任意类型。
键是唯一的，字典只认最后一个赋的键值。

dictionary的方法
D.get(key, 0) #同dict[key]，多了个没有则返回缺省值，0。[]没有则抛异常
D.has_key(key) #有该键返回TRUE，否则FALSE
D.keys() #返回字典键的列表
D.values()
D.items()

D.update(dict2) #增加合并字典
D.popitem() #得到一个pair，并从字典中删除它。已空则抛异常
D.clear() #清空字典，同del dict
D.copy() #拷贝字典
D.cmp(dict1,dict2) #比较字典，(优先级为元素个数、键大小、键值大小)
#第一个大返回1，小返回-1，一样返回0

dictionary的复制
dict1 = dict #别名
dict2=dict.copy() #克隆，即另一个拷贝。
```
### 2. 遍历字典
#### （1）遍历key值
```
>>> a
{'a': '1', 'b': '2', 'c': '3'}
>>> for key in a:
print(key+':'+a[key])
-------------------------------
>>> for key in a.keys():
print(key+':'+a[key])
# 在使用上，for key in a和 for key in a.keys():完全等价。
```
#### （2）遍历value值
```
>>> for value in a.values():
print(value)
```
#### （3）遍历字典项
```
>>> for kv in a.items():
print(kv)

('a', '1')
('b', '2')
('c', '3')
```
#### （4）遍历字典健值
```
>>> for key,value in a.items():
print(key+':'+value)
--------------------------------
>>> for (key,value) in a.items():
print(key+':'+value)
```
***
## 字符串操作
Source:
- [字符串操作](https://blog.csdn.net/zhu_liangwei/article/details/7931701)<br>

```
str = “Hello My friend”
字符串是一个整 体。如果你想直接修改字符串的某一部分，是不可能的。但我们能够读出字符串的某一部分。
子字符串的提取
str[:6]
字符串包含 判断操作符：in，not in
“He” in str
“she” not in str

string模块，还提供了很多方法，如
S.find(substring, [start [,end]]) #可指范围查找子串，返回索引值，否则返回-1
S.rfind(substring,[start [,end]]) #反向查找
S.index(substring,[start [,end]]) #同find，只是找不到产生ValueError异常
S.rindex(substring,[start [,end]])#同上反向查找
S.count(substring,[start [,end]]) #返回找到子串的个数

S.lowercase()
S.capitalize() #首字母大写
S.lower() #转小写
S.upper() #转大写
S.swapcase() #大小写互换

S.split(str, ‘ ‘) #将string转list，以空格切分
S.join(list, ‘ ‘) #将list转string，以空格连接

处理字符串的内置函数
len(str) #串长度
cmp(“my friend”, str) #字符串比较。第一个大，返回1
max(‘abcxyz’) #寻找字符串中最大的字符
min(‘abcxyz’) #寻找字符串中最小的字符

```
***
