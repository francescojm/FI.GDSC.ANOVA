


##### progress bar
my.Pbar.create<-function(max=max){
  pb <- txtProgressBar(min = 0, max = max, style = 3)
  return(pb)
}

my.Pbar.progres<-function(pb,i){
  setTxtProgressBar(pb, i)
}

my.Pbar.close<-function(pb){
  Sys.sleep(0.1)
  close(pb)
}

