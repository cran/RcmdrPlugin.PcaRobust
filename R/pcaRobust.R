#'pcaRobust
#'
#'calculate principal component analysis Robust with ROBPCA method
#'
#'
#'@name pcaRobust
#'@import Rcmdr
#'@import rrcov
#'@import robustbase
#'@import tkrplot
#'@importFrom graphics abline text
#'@importFrom stats window
#'@importFrom utils browseURL
#'@export
#'@examples
#'library(Rcmdr)
#'#pcaRobust()

pcaRobust <- function(){
  tt<-tktoplevel()
  tcltk::tktitle(tt)<-"Robust Principal Component Analysis"

  tl<-tklistbox(tt,height=7,width=20,selectmode="multiple",yscrollcommand=function(...)tkset(scr,...),xscrollcommand=function(...)tkset(scr,...),background="white")
  number<-tklabel(tt,text="Number of components :",fg="blue")
  dataSet<-as.data.frame(eval(parse(text=activeDataSet())))

  # require(foreign)
  list<-colnames(dataSet)
  number<-length(list)
  for (i in (1:length(list)))
  {
    tkinsert(tl,"end",list[i])
  }
  tkselection.set(tl,0)

  pcHubert<-function(){
    Myhscale <- 1.5    # Horizontal scaling
    Myvscale <- 1.5    # Vertical scaling
    alPha=as.numeric(tclvalue(cf)) #confidence interval value
    gettk = as.numeric(tkcurselection(tl))+1
    index<-c()
    for(i in  (1:length(gettk))){
      index[i]<-gettk[i]
    }
    newData<-data.matrix(dataSet[,index])
    pc<-PcaHubert(newData,aplha=alPha,k=length(list),kmax=length(list),mcd=T)
    print(summary(pc))
    cbVal1 <- as.character(tclvalue(cbValue1))
    cbVal2 <- as.character(tclvalue(cbValue2))
    cbVal3 <- as.character(tclvalue(cbValue3))
    cbVal4 <- as.character(tclvalue(cbValue4))
    cbVal5 <- as.character(tclvalue(cbValue5))
    cbVal6 <- as.character(tclvalue(cbValue6))
    cbVal7 <- as.character(tclvalue(cbValue7))

    if (cbVal1=="1"){
      print(pc@loadings)
      }

    if (cbVal2=="1"){
      print(pc@scores)
      }

    if (cbVal3=="1"){
      tt3 <- tktoplevel()
      tkwm.title(tt3,"Screeplot")
      plotFunction <- function()
      {
        screeplot(pc,type="lines",col=3)
      }
      CopyToClip <- function()
      {
        tkrreplot(img)
      }
      img <- tkrplot(tt3,fun=plotFunction,hscale=Myhscale,vscale=Myvscale)
      copy.but <- tcltk::tkbutton(tt3,text="Copy to Clipboard",command=CopyToClip)
      tkgrid(img)
      tkgrid(copy.but)
    }

    if (cbVal4=="1"){
      tt4 <- tktoplevel()
      tkwm.title(tt4,"Biplot")
      plotFunction1 <- function()
      {
        biplot(pc,cex=0.8)
      }
      CopyToClip1 <- function()
      {
        tkrreplot(img1)
      }
      img1 <- tkrplot(tt4,fun=plotFunction1,hscale=Myhscale,vscale=Myvscale)
      copy.but1 <- tcltk::tkbutton(tt4,text="Copy to Clipboard",command=CopyToClip1)
      tkgrid(img1)
      tkgrid(copy.but1)
    }

    if (cbVal5=="1"){
      tt5 <- tktoplevel()
      tkwm.title(tt5,"Orthogonal Distances Plot")
      plotFunction2 <- function()
      {
        y<-pc@od
        x = 1:length(y)
        scaley = (mean(y)-min(y))/5
        xplot = c(x)
        yplot =c(y+scaley)
        plot(xplot,yplot,type="n",xlab = "Number of objects", ylab = "Orthogonal Distance")
        text(x,y+scaley,1:length(y),font=5 , cex=0.6)
        text(x,y,".",font=5, cex=3.6, col=2)
      }
      CopyToClip2 <- function()
      {
        tkrreplot(img2)
      }
      img2 <- tkrplot(tt5,fun=plotFunction2,hscale=Myhscale,vscale=Myvscale)
      copy.but2 <- tcltk::tkbutton(tt5,text="Copy to Clipboard",command=CopyToClip2)
      tkgrid(img2)
      tkgrid(copy.but2)
    }

    if (cbVal6=="1"){
      tt6 <- tktoplevel()
      tkwm.title(tt6,"Score Distances Plot")
      plotFunction3 <- function()
      {
        y<-pc@sd
        x = 1:length(y)
        scaley = (mean(y)-min(y))/5
        xplot = c(x)
        yplot =c(y+scaley)
        plot(xplot,yplot,type="n",xlab = "Number of objects", ylab = "Score Distance")
        text(x,y+scaley,1:length(y),font=5 , cex=0.6)
        text(x,y,".",font=5, cex=3.6, col=2)
      }
      CopyToClip3 <- function()
      {
        tkrreplot(img3)
      }
      img3 <- tkrplot(tt6,fun=plotFunction3,hscale=Myhscale,vscale=Myvscale)
      copy.but3 <- tcltk::tkbutton(tt6,text="Copy to Clipboard",command=CopyToClip3)
      tkgrid(img3)
      tkgrid(copy.but3)
    }

    if (cbVal7=="1"){
      tt7 <- tktoplevel()
      tkwm.title(tt7,"Diagnostic Plot")
      plotFunction4 <- function()
      {
        x<-vector()
        y<-vector()
        h<-vector()
        v<-vector()
        y<-pc@od
        x<-pc@sd
        h<-pc@cutoff.od
        v<-pc@cutoff.sd
        scaley = (mean(y)-min(y))/5
        xplot = c(v,x)
        yplot =c(h,y+scaley)
        plot(xplot,yplot,type="n",xlab = "Score Distance", ylab = "Orthogonal Distance")
        text(x,y+scaley,1:nrow(dataSet),font =5, cex=0.6)
        text(x,y,".",font=5, cex=3.6, col=2)
        abline(h = h,v=v, col = "blue")
      }
      CopyToClip4 <- function()
      {
        tkrreplot(imgplot4)
      }
      imgplot4 <- tkrplot(tt7,fun=plotFunction4,hscale=Myhscale,vscale=Myvscale)
      tkgrid(imgplot4)
      copy.but4 <- tcltk::tkbutton(tt7,text="Copy to Clipboard",command=CopyToClip4)
      tkgrid(copy.but4)
    }
  }

  cancel<-function(){
    tkdestroy(tt)
  }

  reset<-function(){
    tkselection.clear(tl,0,"end")
    cf<-tclVar("0.75")
    entry.cf<-ttkentry(tt,width="10",textvariable=cf)
    cbValueDefault1<-tclVar("0")
    tkconfigure(cb1,variable=cbValueDefault1)
    cbValueDefault2<-tclVar("0")
    tkconfigure(cb2,variable=cbValueDefault2)
    cbValueDefault3<-tclVar("0")
    tkconfigure(cb3,variable=cbValueDefault3)
    cbValueDefault4<-tclVar("0")
    tkconfigure(cb4,variable=cbValueDefault4)
    cbValueDefault5<-tclVar("0")
    tkconfigure(cb5,variable=cbValueDefault5)
    cbValueDefault6<-tclVar("0")
    tkconfigure(cb6,variable=cbValueDefault6)
    cbValueDefault7<-tclVar("0")
    tkconfigure(cb7,variable=cbValueDefault7)
  }

  help<-function(){
    tkgrab.release(window)
    helpIndex <- file.path(system.file("doc",package="RcmdrPlugin.PcaRobust"),"index.html")
    browseURL(helpIndex)
  }

  #gui
  #setting variable
  dataset<-tklabel(tt,text="Data Set : ",fg="blue")
  variabel<-tklabel(tt,text="Variables\n(select two or more)",fg="blue")
  scr <- tkscrollbar(tt, repeatinterval=5,
                     command=function(...)tkyview(tl,...))

  rb1<-tkradiobutton(tt)
  rb2<-tkradiobutton(tt)
  rbValue<-tclVar("default")
  rbValue<-tclVar("user")
  tkconfigure(rb1,variable=rbValue,value="default")
  tkconfigure(rb2,variable=rbValue,value="user")
  default<-tklabel(tt,text="Default")
  user<-tklabel(tt,text="From user, k = ")

  conf<-tklabel(tt,text="Confidence level\n(between 0.5 and 1)",fg="blue")
  cf<-tclVar("0.75")
  entry.cf<-ttkentry(tt,width="10",textvariable=cf)

  plot<-tklabel(tt,text="Plot :",fg="blue")
  cb3<-tkcheckbutton(tt)
  cbValue3<-tclVar("0")
  tkconfigure(cb3,variable=cbValue3)

  cb4<-tkcheckbutton(tt)
  cbValue4<-tclVar("0")
  tkconfigure(cb4,variable=cbValue4)

  cb5<-tkcheckbutton(tt)
  cbValue5<-tclVar("0")
  tkconfigure(cb5,variable=cbValue5)

  cb6<-tkcheckbutton(tt)
  cbValue6<-tclVar("0")
  tkconfigure(cb6,variable=cbValue6)

  cb7<-tkcheckbutton(tt)
  cbValue7<-tclVar("0")
  tkconfigure(cb7,variable=cbValue7)

  options3<-tklabel(tt,text="Screeplot")
  options4<-tklabel(tt,text="Biplot")
  options5<-tklabel(tt,text="Othogonal Distances Plot")
  options6<-tklabel(tt,text="Score Distances Plot")
  options7<-tklabel(tt,text="Diagnostic Plot")

  options<-tklabel(tt,text="Options :",fg="blue")
  cb1<-tkcheckbutton(tt)
  cbValue1<-tclVar("0")
  tkconfigure(cb1,variable=cbValue1)

  cb2<-tkcheckbutton(tt)
  cbValue2<-tclVar("0")
  tkconfigure(cb2,variable=cbValue2)

  options1<-tklabel(tt,text="Factor Loadings")
  options2<-tklabel(tt,text="Scores")

  img1 <- tclVar()
  tclimg1 <- tcltk::tkimage.create("photo", img1, file = system.file("etc", "ok.gif", package="Rcmdr"))
  img2 <- tclVar()
  tclimg2 <- tcltk::tkimage.create("photo", img2, file = system.file("etc", "cancel.gif", package="Rcmdr"))
  img3 <- tclVar()
  tclimg3 <- tcltk::tkimage.create("photo", img3, file = system.file("etc", "help.gif", package="Rcmdr"))
  img4 <- tclVar()
  tclimg4 <- tcltk::tkimage.create("photo", img4, file = system.file("etc", "reset.gif", package="Rcmdr"))

  help.but<-tcltk::tkbutton(tt,text="   Help   ", image=tclimg3, compound="left", command=help,default = "active",borderwidth = 3)
  reset.but<-tcltk::tkbutton(tt,text="  Reset  ", image=tclimg4, compound="left",command=reset, default = "active",borderwidth = 3)
  ok.but<-tcltk::tkbutton(tt,text="  OK   ", image=tclimg1, compound="left",command=pcHubert,default = "active", borderwidth = 3)
  cancel.but<-tcltk::tkbutton(tt,text="  Cancel  ", image=tclimg2, compound="left", command=cancel,default = "active",borderwidth = 3)

  label1<-tklabel(tt,text=" ")
  label2<-tklabel(tt,text="\n")
  label3<-tklabel(tt,text="")
  label4<-tklabel(tt,text="")
  label6<-tklabel(tt,text="")
  label7<-tklabel(tt,text="")

  #setting the position
  tkgrid(label2,column=0,sticky="w")
  tkgrid(variabel,column=1,row=0,sticky="w")
  tkgrid(tl,scr,column=1,row=3,column=1,rowspan=12,sticky="w")
  tkgrid.configure(scr,rowspan=12,column=4,sticky="nsw")
  tkgrid(label7,row=18)
  tkgrid(conf,row=15,column=1,sticky="w")
  tkgrid(entry.cf,row=16,column=1,sticky="w")
  tkgrid(label1,column=9)
  tkgrid(options,column=1,row=19,sticky="w")
  tkgrid(cb1,row=20,column=4,sticky="e")
  tkgrid(cb2,row=21,column=4,sticky="e")
  tkgrid(options1,column=1,row=20,sticky="w")
  tkgrid(options2,column=1,row=21,sticky="w")
  tkgrid(plot,column=10,columnspan=15,row=19,sticky="w")
  tkgrid(cb3,row=20,column=12,sticky="w")
  tkgrid(cb4,row=21,column=12,sticky="w")
  tkgrid(cb5,row=22,column=12,sticky="w")
  tkgrid(cb6,row=23,column=12,sticky="w")
  tkgrid(cb7,row=24,column=12,sticky="w")
  tkgrid(options3,column=10,row=20,sticky="w")
  tkgrid(options4,column=10,row=21,sticky="w")
  tkgrid(options5,column=10,row=22,sticky="w")
  tkgrid(options6,column=10,row=23,sticky="w")
  tkgrid(options7,column=10,row=24,sticky="w")
  tkgrid(label6,column=17)
  tkgrid(label3,row=25)
  tkgrid(help.but,row=26,column=1,columnspan=15,sticky="w")
  tkgrid(ok.but,row=26,column=1,columnspan=5,sticky="se")
  tkgrid(reset.but,row=26,column=10,columnspan=20,sticky="nsw")
  tkgrid(cancel.but,column=1,row=26,columnspan=15,sticky="e")
  tkgrid(label4,row=27)

}
