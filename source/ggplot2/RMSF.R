library("ggplot2")
#================================================================#
#====   Step 1: �ļ�Ԥ����   ====================================#
#================================================================#
rmsf1 = read.table("D:/Github/DouBeeTwT.github.io/source/ggplot2/rmsf_md1.out") #���ļ�����rmsf1
rmsf2 = read.table("D:/Github/DouBeeTwT.github.io/source/ggplot2/rmsf_md2.out") #���ļ�����rmsf2

rmsf = rbind(rmsf1,rmsf2)                       #�ϲ�rmsf1��rmsf2
rmsf = cbind(rmsf,rep(c("md1","md2"),each=282)) #������Դ��ʶ��
colnames(rmsf) = c("Res","Num","Group")         #����������


#================================================================#
#====   Step 2: �����ʽ����   ==================================#
#================================================================#
t = theme_bw() +                                                    #��ո�ʽ
  theme(panel.grid =element_line(color="grey",linetype="dotted")) + #�޸�������Ϊ��ɫ����
  theme(panel.border = element_blank()) +                           #ɾ�����߿�
  theme(axis.line = element_line()) +                               #����X���Y��
  theme(legend.position=c(0.9,0.9)) +                               #��ǩλ��
  theme(legend.background = element_blank()) +                      #ɾ����ǩ����ɫ
  theme(legend.title=element_blank()) +                             #ɾ����ǩ����
  theme(legend.text = element_text(size = 13)) +                    #�޸ı�ǩ�����С
  theme(axis.title = element_text(size = 15)) +                     #�޸���������������С
  theme(axis.text = element_text(size = 13))                        #�޸�������̶������С


#================================================================#
#====   Step 3: ע�Ͳ�������   ==================================#
#================================================================#
a1 = annotate("rect",xmin=5,xmax=12,ymin=0,ymax=4,fill="black",alpha=0.1)     #���Ӹ�������1
a2 = annotate("rect",xmin=173,xmax=180,ymin=0,ymax=4,fill="orange",alpha=0.1) #���Ӹ�������2
a3 = annotate("rect",xmin=235,xmax=248,ymin=0,ymax=4,fill="purple",alpha=0.1) #���Ӹ�������3
a4 = annotate("rect",xmin=0,xmax=300,ymin=0,ymax=2,fill="white",alpha=0.8)    #������������

line = data.frame(x=c(0,300),y=c(2,2))           #�߶����˵����ݳɱ�
l = geom_hline(yintercept = 2,linetype="dashed") #��������


#================================================================#
#====   Step 4: ��������ͼ   ====================================#
#================================================================#
p = ggplot(rmsf,aes(Res,Num,color=Group)) + geom_line(size = 0.8) + #��������������ͼ
  xlim(0,300) + ylim(0,4) +                                         #X����Y����ʾ��Χ
  xlab("Residues Number") + ylab("RMSF(\uc5)")                      #X����Y�������������
p = p + t + a1 + a2 +a3 + a4 + l                                    #���������ʽt,ע��a1-a4,������l
p                                                                   #��ʾ