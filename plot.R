# Set Working Directory
setwd(dir = "~/Kiryu/201411Mixture/output/")

# animation library
library(animation)

# original data and results of my program
data <- read.table(file="../input/sample_mix.txt", head = T)
results <- read.table(file="../output/results.txt", head = T)

# EMの繰り返し回数がi回の時のデータと混合ガウス分布の円を描く関数
plot_step <- function(i){
  # まず，円を描く関数を定義
  circle <- function(x, y, r){
    th <- seq(0, 2*pi, length.out=100)
    par(new=T)
    plot(x+r*cos(th), y+r*sin(th),
         xlim=c(0,1), ylim=c(0,1),
         type="l", ann=F, asp=1)
  }
  # 最後に描く円のみ，ラベルを印字する
  circle_LAST <- function(x, y, r){
    th <- seq(0, 2*pi, length.out=100)
    par(new=T)
    plot(x+r*cos(th), y+r*sin(th),
         xlim=c(0,1), ylim=c(0,1),
         xlab="x", ylab="y",
         type="l", ann=T, asp=1)
  }
  
  # dataをplotする
  par(new = F)
  plot(data, xlim=c(0,1), ylim=c(0,1), ann = F, asp=1)

  # 円をplotする
  par(new = T)
  k_max <- ((length(results[0,]) - 1)/4) - 1
  for(k in 0:(k_max - 1)){
    circle(results[i,4 * k + 4], 
           results[i,4 * k + 5],
           results[i,4 * k + 3])
  }
  circle_LAST(results[i,4 * k_max + 4], 
              results[i,4 * k_max + 5],
              results[i,4 * k_max + 3])
  par(new = F)
}

animation::saveGIF(
  {for(i in 1:10){plot_step(i)}},
  imgdir = "output",
  movie.name = "results.gif",
  interval=0.1,　ani.width = 600, ani.height = 600
)
