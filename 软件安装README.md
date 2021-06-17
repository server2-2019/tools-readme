# tools-readme
Record readme of often used software
2019.01.11
Foxit PDF Editor 安装, 目的pdf去保护
① 太平洋下载FoxitPDFEditor21_2.2.0.205.zip；
② 解压缩 “FoxitPDFEditor21_2.2.0.205.zip”
③ 双击打开应用程序程序
④ 打开目标pdf => 文档 => 导出页面 => 选择目的地址，指定导出页页面范围

2019.01.17
everything 下载安装
① download Everything-1.4.1.895.x86-Setup from https://www.sogou.com/link?url=hedJjaC291PB9HkP95Hbb1odZoieaM0q_TVVI34IoCfTqOp6wwSpBg..
② 将Everything-1.4.1.895.x86-Setup文件copy到D:\moonbio\software；
③ 双击该安装文件 Everything-1.4.1.895.x86-Setup 进行软件安装


2019.01.22
chrome installation
① download “71.0.3544.2_chrome_installer_win64.exe” from https://repo.fdzh.org/chrome/exe/

2019.02.02
shadowsocks vnp 搭建
① 下载 7z解压软件 “https://7-zip.org/a/7z1806.exe” from “https://sourceforge.net/p/sevenzip/discussion/45797/thread/f1a142fa03/”
  下载后双击安装
② 下载 “ShadowsocksR-4.6.1-win.7z” 并用 7z解压, 下载“OmegaOptions.bak”
③ 下载 “SwitchyOmega_Chromium.crx” from https://github.com/FelisCatus/SwitchyOmega/releases，并将“SwitchyOmega_Chromium.crx”加入chrome扩展文件（工具->扩展程序->将“SwitchyOmega_Chromium.crx”拖入打开的扩展程序页面）
④ 

2019.02.11
enter root 
sudo passwd ruijuan 72019

2019.02.12
FileZilla
--------------------------------------------------------------------------------
# useful links
FileZilla连接ftp服务器(centos)传文件: https://www.jianshu.com/p/b0c146f4ab87
 
--------------------------------------------------------------------------------

① download "FileZilla_3.40.0_win64-setup_bundled" from https://filezilla-project.org/download.php?type=client
② 服务器安装vsftpd
  $ rpm -qa | grep vsftpd # 查是否安装了vsftpd
  $ sudo yum  -y install vsftpd
  $ chkconfig vsftpd on # 设置开机启动
  $ 注释掉/etc/vsftpd/ftpusers文件中的root，因为这一行禁止root登录，然后启动vsftpd服务
    systemctl start vsftpd
  $ 主机：sftp://192.168.10.74 
    端口: 22
    # 注意的是sftp协议默认使用的是22端口，这个协议依赖于sshd服务，所以sshd服务需要启动

  2019.02.20
  google访问助手安装
  ① 下载 360.zip from http://www.ggfwzs.com/ ,解压# 包含google访问插件，及安装说明
  ② 下载 360se10.0.1634.0 from https://browser.360.cn/
  ③ 安装 360se10.0.1634.0
  ④ 安装google访问助手插件： 打开360安全浏览器->右上角管理->管理->将第一步下载的“谷歌访问助手_v2.3.0.crx”拖入打开的管理页面
  ⑤ 设置永久有效：右上角打开菜单->设置->基本设置->修改主页为“https://2018.hao245.com/”
  ⑥ 关闭 360浏览器，从新打开，
  ⑦ 可用的google scholar 如 https://scholar.google.co.jp/ # 香港连接不可用

  2019.02.22
  MobaXterm更新
  下载 portabl edition “MobaXterm_Portable_v11.1.zip” from https://mobaxterm.mobatek.net/download-home-edition.html
  解压 MobaXterm_Portable_v11.1.zip 

  2019.03.04
  MEGA X 
  download MEGA_X_10.0.5_win64_setup.exe from /work/workspace/yuanjie/download/MEGA_X_10.0.5_win64_setup.exe to D:\moonbio\software
  Destination location: C:\Program Files\MEGA-X

  2019.03.13
  添加Feedly Subscribe Button 到chrome using the URL as following:
  https://chrome.google.com/webstore/detail/feedly-subscribe-button/gbbnddjfcllebfcnihfgmdplgaiejepc

  XMind 
  XMind download from https://www.xmind.net/download/

  Feed43 创建RSS
  https://feed43.com/feed.html?name=5225417662588288

  Mendeley
  download Mendeley-Desktop-1.19.3-win32 from https://www.mendeley.com/download-desktop/#downloading
  add Mendeley Importer to chrome: https://chrome.google.com/webstore/detail/mendeley-importer/dagcmkpagjlhakfdhnbomgmjdpkdklff
  TIPS:
  1.同步设置，由于每个人只有2G的容量使用空间，所以这里不同步pdf文件，只同步相应的设置信息和bibtex信息。 
就是不要在同步attached files哪儿打钩就可以了。 
  2.文件管理，这儿文件管理我有自己的一套方案，因此只需要设置文献的存储路径就好了。其他两项千万别用，不然软件会自动修改文件的储存路径和文件的命名方式，最主要的还TM不智能。对于有强迫症的我就只设置第一项了。
  3.bibtex的设置，这个设置很重要，我们对pdf文献进行的一些关联信息的修改都存储这bibtex文件之中。但是这儿的设置很简单，我们只需要设置文件的存储方式和文件的位置就好了。我的设置如下 
  注意：初次在新电脑上进行同步工作的时候，首先进行上述配置，然后不要同步，而是手动导入bib文件，在File->import->BibTex，然后导入原来的bib文件，这样同步工作就ok了。

  install Sci-Hub extension to chrome
  To do this:
  Download the extension and unpack it. You get the "Sci-Hub" folder with code.
  Open Chrome and navigate to chrome://extensions, or just open the menu -> settings -> extensions.
  Check the developer mode in upper right.
  If you have old version of Sci-Hub extension installed, remove it
  Click "Load unpacked extension" button
  Highlight the folder "Sci-Hub" (do not enter it) and click "Open".
  Done. Go back to sci-hub.tw and use search!

 2019.03.26
 破解版XMIND 安装
 从雅俊电脑copy “Xmind Zen” 安装软件到 D:\moonbio\software 双击安装 + 破解补丁

 2019.04.04
 BaiduNetdisk_6.7.2.16 百度网盘客户端安装
 下载BaiduNetdisk_6.7.2.16 from 百度网盘，存放地址 D:\moonbio\software，双击安装
 
 安装youtube downloader
 安装说明参考：https://www.wikihow.com/Download-YouTube-Videos-in-Chrome
 1. download "youtube_video_downloader_16.1.2_chrome.zip" from https://addoncrop.com/youtube_video_downloader/ to D:\moonbio\software\chromeYoutubeDownloader;
 2. double click "youtube_video_downloader_16.1.2_chrome.zip" 解压到 “youtube_video_downloader_16.1.2_chrome” 中
 2. 打开chrome extension 页面，点击“加载已解压的扩展程序”， 加载插件

 印象笔记chrome剪藏插件安装
 下载 Evernote-Web-Clipper_v7.0.2.crx 到 D:\moonbio\software\印象笔记 from https://www.yinxiang.com/webclipper/install/
 添加下载的插件到chrome扩展程序

 添加chrome插件to-do-list
 1. open https://chrome.google.com/webstore/detail/todoist-to-do-list-and-ta/jldhpllghnbhlbpcmnajkpdmadaolakh?hl=zh-CN中添加；
 2. click “添加至chrome”

 Scrivener文章编写软件
 1. download Scrivener_16734.zip from http://www.zdfans.com/html/16734.html
 2. 双击解压 Scrivener_16734.zip，
 3. 破解：安装完成后到桌面找到软件，单击鼠标右键——打开文件位置(C:\ProgramData\Microsoft\Windows\Start Menu\Programs\Scrivener)=> 将数据包中fix文件夹下提供的破解补丁复制到上一步打开的文件目录下替换掉原文件即可破解

 github desktop install
 1. download "GitHubDesktopSetup" from https://desktop.github.com/ to D:\moonbio\software\github
 2. double click "GitHubDesktopSetup"

 2019.05.07
 zhangboVPN设置 (Shandowsocks cordcloud)

 1. 注册（CordCloud）
    ① 打开：https://www.cordcloud.cc/auth/register1?code=BBC4ETGpVYfFzCEmtoNQI11ijqLaJakR；
    ② 姓： moon; 名：bio
      email: moonb8467@gmail.com (邮箱信息 pwd：moonb8467bio；生日：1997.07.01)
      pwd: moonb8467cordcloud
    ③ 左侧导航栏 充值 => 充值150元；
    ④ 左侧导航栏 商店 => 流量套餐1750G（可叠加）150.0元
    ⑤ 左侧导航栏 首页 => SHADOWSOCKSR, WINDOWS => 下载 《ssr-win.7z》 到 D:\moonbio\software\DrZhangVPN => 解压 《ssr-win.7z》 到 ssr-win
    ⑥ open D:\moonbio\software\DrZhangVPN\ssr-win => 双击 ShadowsocksR 打开Shadowsocks
    ⑦ Shadowsocks 设置： (推荐)右键小飞机--服务器--SSR服务器订阅设置，将订阅地址设置为下面的地址，其他参数留空，确定之后再更新 SSR 服务器订阅。
                         然后选择一个合适的服务器(选择时请右击小飞机选择)，务必关闭负载均衡，系统代理模式选”全局模式”
      备注： 如出现无法更新订阅情况请将www.cordcloud.cc替换为www.cordcloud.online
            普通端口地址：https://www.cordcloud.cc/link/FObAFQlJU9MwOYpJ?mu=0
            443端口地址：https://www.cordcloud.cc/link/FObAFQlJU9MwOYpJ?mu=443
            80端口地址：https://www.cordcloud.cc/link/FObAFQlJU9MwOYpJ?mu=80
            543端口地址：https://www.cordcloud.cc/link/FObAFQlJU9MwOYpJ?mu=543
      1）Mode： global
      2) Proxy rule： bypass LAN & China
      3）Servers: 选择合适的 SSR node
      4）Servers Subscribe: 
        （1）Subscribe setting 添加
            普通端口地址：https://www.cordcloud.cc/link/FObAFQlJU9MwOYpJ?mu=0
            443端口地址：https://www.cordcloud.cc/link/FObAFQlJU9MwOYpJ?mu=443
            80端口地址：https://www.cordcloud.cc/link/FObAFQlJU9MwOYpJ?mu=80
            543端口地址：https://www.cordcloud.cc/link/FObAFQlJU9MwOYpJ?mu=543
        （2）Update Subscribe SSR node

  2020.02.12更新
  下载，解压，运行程序(建议使用管理员方式启动)，您有三种方式导入所有节点
  (1)下载这个（普通端口）或者这个（单端口多用户），右键小飞机 服务器 -- 从配置文件导入服务器，选择这个文件，
  (2)点击这里（普通端口）或者这个(单端口多用户），然后右键小飞机 -- 从剪贴板复制地址
  (3)(推荐)右键小飞机--服务器--SSR服务器订阅设置，将订阅地址设置为下面的地址，其他参数留空，确定之后再更新 SSR 服务器订阅。
  然后选择一个合适的服务器(选择时请右击小飞机选择)，务必关闭负载均衡，系统代理模式选”全局模式”，然后即可上网。
  
  SSR 订阅地址(推荐)：
  如出现无法更新订阅情况请将www.cordcloud.club替换为www.cordcord.xyz
  普通端口地址：https://www.cordcord.xyz/link/FObAFQlJU9MwOYpJ?mu=0
  443端口地址：https://www.cordcord.xyz/link/FObAFQlJU9MwOYpJ?mu=443
  80端口地址：https://www.cordcord.xyz/link/FObAFQlJU9MwOYpJ?mu=80
  543端口地址：https://www.cordcord.xyz/link/FObAFQlJU9MwOYpJ?mu=543

  Clash
  https://www.cordcord.xyz/link/SdVwtMGhTOr0ij1s?mu=4
  https://tgbot.lbyczf.com/surge2clash?url=https://www.cordcord.xyz/link/SdVwtMGhTOr0ij1s?mu=4

  SSD


  snlVNP
  email:ruijuan6(ruijuan6(3))
  1. 注册账户（https://www.nsl-net.com/auth/register?code=mLo0）
     ruijuan6（ruijuan6（3））
     name:bio
  2. 充值： 进入账户左侧列表充值30元；
  3. 套餐购买：进入账户左侧列表“套餐购买”，vip3 30元标准不限制流量
  4. 软件下载设置：
    （1）左侧列表“用户中心”=> SSR+WINDOWS 下载客户端《ssr-win.7z》 到 D:\moonbio\software\nslVPN；
    （2）加压《ssr-win.7z》 => 双击“ShadowsocksR”
     (3) 右键小飞机 => Server Subscribe -> Server Setting, 添加：
         https://www.nsl-net.com/link/EO26gZVtolQskN8N?mu=0
         https://www.nsl-net.com/link/EO26gZVtolQskN8N?mu=1
         右键小飞机 => Server Subscribe -> Update Subscribe SSR node,
    （4）Mode PAC, 然后选择一个合适的服务器，代理规则选“绕过局域网和大陆”，然后即可上网
         右键小飞机 => Server Subscribe -> Update Subscribe SSR node,
  android NSL 安装
  1. 下载 “ssr-android” 到 D:\moonbio\software\android-nsl
  2. 打开App，点击右下角的add号图标
  3. 添加/升级 SSR订阅
  4. 添加订阅地址，输入下方订阅地址后确定
  5. 订阅出现系统自带的与NSL Network，请把系统自带的无效订阅左滑删除（自带影响订阅更新速度）
  6. 点击确定并升级
  7. 点击选择任意节点， 路由选择：略过区域网路以及中国大陆
  8. 点击右上角的纸飞机图标即可连接
  9. 备用导入节点方法：在手机上默认浏览器中点击普通端口链接或者单端口多用户链接，然后点击确定
     普通节点： https://www.nsl-net.com/link/EO26gZVtolQskN8N?mu=0
     单端口节点： https://www.nsl-net.com/link/EO26gZVtolQskN8N?mu=1

   graphviz 安装
   ref: https://feedly.com/i/entry/tpNXTg4zkUnw5NWQrD2jbXoF6soRjpnPjWqAP2Z37Io=_16a9053a503:10d39e1:d10d5a25
   ref: https://graphviz.gitlab.io/_pages/Download/Download_source.html
   1. download graphviz.tar.gz from https://graphviz.gitlab.io/_pages/Download/Download_source.html to /mnt/d/home/ruijuan/workflows/software
   2. cd /mnt/d/home/ruijuan/workflows/software
      tar zxvf graphviz.tar.gz
      cd graphviz-2.40.1/
      ./configure
      make 
      make install 
      

  2019.05.20
  NextCloud 安装
  1. download “Nextcloud-2.5.2-setup” from https://nextcloud.com/install/#install-clients to D:\moonbio\software\nextCloud
  2. 双击 “Nextcloud-2.5.2-setup” 进行安装

  2019.05.22
  FSCapture 截图工具
  1. download "FSCaptureSetup90.exe" from https://www.faststone.org/FSCaptureDownload.htm to D:\moonbio\software
  2. double click "FSCaptureSetup90.exe"
  
  2019.06.03
  nextcloud reset password for the account
  admin account: x
  users -> change password -> click "Enter"

  2019.06.11
  ref: https://github.com/gildas-lormeau/SingleFileZ
  SingleFileZ – 网摘新工具：打包压缩完整网页
  SingleFileZ添加到chrome扩展
  1. open Chrome web store: https://chrome.google.com/webstore/detail/singlefilez/offkdfbbigofcgdokjemgjpdockaafjg
  2. 添加扩展到chrome
  3. chrome -> 更多工具 -> 扩展程序 -> 详细信息 -> 打开“允许访问文件网址”
  Chrome: Install SingleFileZ and enable the option "Allow access to file URLs" in the details page of the extension (chrome://extensions/?id=offkdfbbigofcgdokjemgjpdockaafjg). Otherwise, start the browser with the switch "--allow-file-access-from-files".

  zenflowchart 在线流程图制作软件
  ref：https://www.zenflowchart.com/?ref=appinn

  2019.Jun17-09:03
  nextcloud restart command:
  nohup python /mnt/d/linux/M/www/Django/moonbio/manage.py runserver 0:9000 &
  domain name: nextcloud.mn

  医家汇
  https://www.ediecogroup.com/method_article_detail/76/
  zrjlyq@126.com
  8yjh

  2019.Jul29-15:44
  新服务器账号密码
  rock.cluster.mn
  ssh 端口 22
  朱瑞娟：zhurj 
  passwd: moonbio1

  192.168.2.74
  ruijuan
  (7)2019

Moonbio Tech
moon@123

# csdn
13632203310@pwd8+2019

2019.Aug15-09:48
worktile账户密码
zhurj@moonbio.com+pwd8

2019.Aug27-08:46
URL: https://hub.docker.com/
Docker
shitou62019
ruijuan6@gmail.com
8docker

2019.Sep23-13:19
chrome json handle 插件安装
1. 访问http://jsonhandle.sinaapp.com/下载 JSON-Handle
2. chrome 设置 -> 更多工具 -> 扩展程序 -> JSON-handle_0.6.1.crx 拖到扩展程序中

2019.Oct11-10:25
galaxy
1. 账户 ： 126
2. 8galaxy

2019.Oct29-10:10
github
server2-2019
847052660@qq.com
8server2-2019


http://moondisk.mn:5000/
4ppt
ppt@moon

\\moondisk.mn
zhurj
Hello123

2020.Jan13-16:11
reddit.com
ruijuan6@gmail.com
stone6_2020
8

ruijuan6@gmail.com
zqzrj[a-z]{3}\d{4}

https://www.oracle.com/index.html
c:126
p:8

2020.Apr03-14:51
PUBMLST
https://pubmlst.org/bigsdb
USERNAME: Stone6
passwd: 8pubmlst
email: ruijuan6@gmail.com

2020.06.16
删除环境变量
C:\Users\MoonBiotech\AppData\Local\Microsoft\WindowsApps

2020.06.30 添加环境变量
windows桌面搜索框中输入：“系统” => 点击“系统”（蓝色小电脑，非系统信息）=> 高级系统设置 => 高级 => 环境变量 => 系统环境变量 => Path => 新建  输入 需要添加的环境变量

2020.07.03
Feedly account
126.com
8

supermoon2019@outlook.com

2020.08.17
向日葵远程控制
https://pc.qq.com/detail/8/detail_13068.html
购买了入门版本： 截止日期2021.08.17


2020.09.08
修改zhouyj密码
moon2020

2020.09.15
坚果云
126.com
8jianguoyun
stone6

2020.09.25
https://forum.biobakery.org/u/account-created
gmail.com
stone
joan
8humann
zqxxxwhl1024

有比例的venn图
参考：https://omics.pnl.gov/software/venn-diagram-plotter
网页作图：http://www.biovenn.nl
软件作图：Venn Diagram Plotter
D:\moonbio\software\VennDiagramPlotter_Installer_0


ppt 模板
timeline info
https://www.free-powerpoint-templates-design.com/free-powerpoint-diagrams-design/free-powerpoint-timeline-diagrams/page/2/


2021.03.15 购买processon 个人版 159元


