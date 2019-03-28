library("ggplot2")
#================================================================#
#====   Step 1: 文件预处理   ====================================#
#================================================================#
rmsf1 = read.table("D:/Github/DouBeeTwT.github.io/source/ggplot2/rmsf_md1.out") #从文件读入rmsf1
rmsf2 = read.table("D:/Github/DouBeeTwT.github.io/source/ggplot2/rmsf_md2.out") #从文件读入rmsf2

rmsf = rbind(rmsf1,rmsf2)                       #合并rmsf1和rmsf2
rmsf = cbind(rmsf,rep(c("md1","md2"),each=282)) #添加来源标识列
colnames(rmsf) = c("Res","Num","Group")         #重命名列名


#================================================================#
#====   Step 2: 主题格式设置   ==================================#
#================================================================#
t = theme_bw() +                                                    #清空格式
  theme(panel.grid =element_line(color="grey",linetype="dotted")) + #修改网格线为灰色虚线
  theme(panel.border = element_blank()) +                           #删除外层边框
  theme(axis.line = element_line()) +                               #添加X轴和Y轴
  theme(legend.position=c(0.9,0.9)) +                               #标签位置
  theme(legend.background = element_blank()) +                      #删除标签背景色
  theme(legend.title=element_blank()) +                             #删除标签标题
  theme(legend.text = element_text(size = 13)) +                    #修改标签字体大小
  theme(axis.title = element_text(size = 15)) +                     #修改坐标轴标题字体大小
  theme(axis.text = element_text(size = 13))                        #修改坐标轴刻度字体大小


#================================================================#
#====   Step 3: 注释部件构建   ==================================#
#================================================================#
a1 = annotate("rect",xmin=5,xmax=12,ymin=0,ymax=4,fill="black",alpha=0.1)     #添加高亮区域1
a2 = annotate("rect",xmin=173,xmax=180,ymin=0,ymax=4,fill="orange",alpha=0.1) #添加高亮区域2
a3 = annotate("rect",xmin=235,xmax=248,ymin=0,ymax=4,fill="purple",alpha=0.1) #添加高亮区域3
a4 = annotate("rect",xmin=0,xmax=300,ymin=0,ymax=2,fill="white",alpha=0.8)    #弱化无用区域

line = data.frame(x=c(0,300),y=c(2,2))           #线段两端点数据成表
l = geom_hline(yintercept = 2,linetype="dashed") #添加虚线


#================================================================#
#====   Step 4: 整合与作图   ====================================#
#================================================================#
p = ggplot(rmsf,aes(Res,Num,color=Group)) + geom_line(size = 0.8) + #现有数据作线型图
  xlim(0,300) + ylim(0,4) +                                         #X轴与Y轴显示范围
  xlab("Residues Number") + ylab("RMSF(\uc5)")                      #X轴与Y轴坐标标题内容
p = p + t + a1 + a2 +a3 + a4 + l                                    #添加主题格式t,注释a1-a4,辅助线l
p                                                                   #显示