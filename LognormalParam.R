
#Distribución lognormal para desviación estándar: 

#SD<-vector de la desviación estándar medida

#sel_bin_vel: número para seleccionar las desviaciones estándar por bin de velocidad

#method: Puede ser Mom (Momentos) ó MLE (Máxima Verosimilitud) 

#Vel: vector de velocidad asociado a la desviación estándar

#ADVERTENCIA: MLE requiere más esfuerzo computacional (Con sel_bin_vel activado, va mucho más rápido!)

LognormalParam<-function(SD,Vel, sel_bin_vel=NULL, method){
  
  Vel_bin<-floor(Vel)
  
  require("ggplot2"); require("ggthemes")
  
  if(is.vector(SD)==FALSE | is.vector(Vel)==FALSE){
    stop("SD y Vel deben ser vectores")
  }
  
  if(any(SD[is.na(SD)==FALSE]>30)==TRUE | any(Vel[is.na(Vel)==FALSE]>90)){
    stop("SD y Vel deben ser filtrados. Valores superiores a 30 (90) encontrados")
  }
  
  if(is.numeric(SD)==FALSE | is.numeric(Vel)==FALSE){
    stop("SD y Vel deben ser numéricos")
  }
  
  if(is.character(method)==FALSE){
    stop("method es un argumento inválido")
  }
  
  
  
  if(method=="Mom"){
    
    SD_df<-as.data.frame(matrix(nrow = length(SD), ncol = 2)); SD_df[,1]<-SD; SD_df[,2]<-Vel_bin
    
    SD_df<-SD_df[is.na(SD_df$V1)==FALSE,]
    
    SD[SD=="NaN"]<-NA
    
    SD<-SD[is.na(SD)==FALSE]
    
    SD1<-round(SD,1)
    
    
    
    if(is.numeric(sel_bin_vel)==TRUE){
      SD<-SD_df$V1[SD_df$V2==sel_bin_vel]
      SD1<-round(SD,1)
    }
    else{
      SD<-SD
      SD1<-SD1
    }
    
    suma1<-sum(SD); suma2<-sum(SD^2)
    
    n<-length(SD)
    
    mu<-(-log(suma2)/2)+2*log(suma1)-((3/2)*log(n))
    
    s<-log(suma2)-2*log(suma1)+log(n)
    
    sigma<-sqrt(s)
    
    
    
    #Cálculo de frecuencias:
    
    SD2<-as.vector(table(SD1))
    suma3<-sum(SD2)
    frecuencia<-SD2/suma3
    
    SDbinmedio<-as.vector(tapply(SD,SD1, mean))
    
    frecuencia1<-vector(mode = "numeric", length = length(frecuencia))
    
    for(i in 1:length(frecuencia)){
      frecuencia1[i]<-exp(-((log(SDbinmedio[i])-mu)^2/(2*sigma^2)))*(1/(sqrt(2*pi)*sigma*SDbinmedio[i]))
    }
    
    frecuencia1<-frecuencia1/sum(frecuencia1)
    
    suma_chi<-0
    
    for(i in 1:length(frecuencia1)){
      suma_chi<-suma_chi+(frecuencia[i]-frecuencia1[i])^2
    }
    
    media<-mean(frecuencia, na.rm = TRUE)
    suma11<-0
    
    for(i in 1:length(frecuencia1)){
      suma11<-suma11+(frecuencia[i]-media)^2
    }
    
    R<-1-(suma_chi/suma11)
    
    Chi_Cuadrado<-suma_chi/(length(SDbinmedio)-2)
    
    RMSE<-(suma_chi/length(SDbinmedio))^0.5
    
    Coeficientes<-list(mu=mu, sigma=sigma, Chi_Cuadrado=Chi_Cuadrado, RMSE=RMSE, R_2=R)
    
    z<-as.numeric(names(tapply(SD,SD1,mean)))
    
    #Gráfico de la desviación estándar: (el gráfico coincide! pero el resultado NO!!)
    
    plot(x=SDbinmedio, y=frecuencia, col = "red", main = "Distribución Log-Normal de SD", xlab = "Desviación estándar (m/s)", ylab = "Frecuencia", pch=2, lwd=2, ylim = c(0,(max(frecuencia1)+0.05)))
    legend("topright", col = c("blue", "red"), legend = c("Ajuste log-normal", "Datos medidos"), lwd = 2, bty = "n", border = "black", box.lwd = 2, box.col = "black", box.lty = "l")
    lines(x=SDbinmedio, y=frecuencia1, col="blue", lwd=2)
    
    
    df<-as.data.frame(matrix(nrow = length(frecuencia), ncol = 3))
    
    df[,1]<-SDbinmedio; df[,2]<-frecuencia; df[,3]<-frecuencia1
    
    p<-ggplot(data = df, aes(x=V1))+
      geom_point(aes(x=V1,y=V2),size=4, shape=21, fill="blue")+
      geom_line(aes(y=V3, colour="blue"), lwd=1.5)+
      scale_x_continuous(name = "Desviación Estándar (m/s)") +
      scale_y_continuous(name = "Frecuencia") +
      ggtitle(paste("Distribución Log-Normal para Vbin=", sel_bin_vel, "(m/s) por Momentos")) +
      theme_economist() +
      theme(legend.position = c(0.75,0.75),
            legend.background = element_rect(fill = "white",linetype = "solid", colour="lightblue", size = 2),
            plot.title = element_text(size = 18, margin = margin(10,10,10,10)),
            text = element_text(),
            axis.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 9, face = "bold"),
            legend.title=element_blank())+
      scale_colour_discrete(labels=c("Ajuste log-normal"))
    
    
    Coeficientes<-list(mu=mu, sigma=sigma, Chi_Cuadrado=Chi_Cuadrado, RMSE=RMSE, R_2=R, graph=p)
    
    
    return(Coeficientes)
    
  }
  
  
  if(method=="MLE"){
    
    SD_df<-as.data.frame(matrix(nrow = length(SD), ncol = 2)); SD_df[,1]<-SD; SD_df[,2]<-Vel_bin
    
    SD_df<-SD_df[is.na(SD_df$V1)==FALSE,]
    
    SD[SD=="NaN"]<-NA
    
    SD<-SD[is.na(SD)==FALSE]
    
    SD1<-round(SD,1)
    
    if(is.numeric(sel_bin_vel)==TRUE){
      SD<-SD_df$V1[SD_df$V2==sel_bin_vel]
      SD1<-round(SD,1)
    }
    else{
      SD<-SD
      SD1<-SD1
    }
    
    suma1<-0
    
    n<-length(SD)
    
    for(i in 1:n){
      suma1<-suma1+log(SD[i])
    }
    
    mu<-suma1/n
    
    suma2<-vector(mode = "numeric", length = n)
    suma3<-vector(mode = "numeric", length = n)
    
    #REVISAR ESTE CÁLCULO AQUÍ! HAY FALLO!
    
    for(i in 1:n){
      for(j in i:n){
        suma3[i]<-suma3[i]+log(SD[j])
      }
      suma2[i]<-suma2[i]+(log(SD[i])-(suma3[i]/n))^2
    }
    
    suma2<-sum(suma2)
    
    sigma<-sqrt(suma2/n)
    
    #Cálculo de frecuencias:
    
    SD2<-as.vector(table(SD1))
    suma3<-sum(SD2)
    frecuencia<-SD2/suma3
    
    SDbinmedio<-as.vector(tapply(SD,SD1, mean))
    
    frecuencia1<-vector(mode = "numeric", length = length(frecuencia))
    
    for(i in 1:length(frecuencia)){
      frecuencia1[i]<-exp(-((log(SDbinmedio[i])-mu)^2/(2*sigma^2)))*(1/(sqrt(2*pi)*sigma*SDbinmedio[i]))
    }
    
    frecuencia1<-frecuencia1/sum(frecuencia1)
    
    suma_chi<-0
    
    for(i in 1:length(frecuencia1)){
      suma_chi<-suma_chi+(frecuencia[i]-frecuencia1[i])^2
    }
    
    media<-mean(frecuencia, na.rm = TRUE)
    suma11<-0
    
    for(i in 1:length(frecuencia1)){
      suma11<-suma11+(frecuencia[i]-media)^2
    }
    
    R<-1-(suma_chi/suma11)
    
    Chi_Cuadrado<-suma_chi/(length(SDbinmedio)-2)
    
    RMSE<-(suma_chi/length(SDbinmedio))^0.5
    
    z<-as.numeric(names(tapply(SD,SD1,mean)))
    
    df<-as.data.frame(matrix(nrow = length(frecuencia), ncol = 3))
    
    df[,1]<-SDbinmedio; df[,2]<-frecuencia; df[,3]<-frecuencia1
    
    p<-ggplot(data = df, aes(x=V1))+
      geom_point(aes(x=V1,y=V2),size=4, shape=21, fill="blue")+
      geom_line(aes(y=V3, colour="blue"), lwd=1.5)+
      scale_x_continuous(name = "Desviación estándar (m/s)") +
      scale_y_continuous(name = "Frecuencia") +
      ggtitle(paste("Distribución Log-Normal para Vbin=", sel_bin_vel, "(m/s) por MLE")) +
      theme_economist() +
      theme(legend.position = c(0.75,0.75),
            legend.background = element_rect(fill = "white",linetype = "solid", colour="lightblue", size = 2),
            plot.title = element_text(size = 18, margin = margin(10,10,10,10)),
            text = element_text(),
            axis.title = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 9, face = "bold"),
            legend.title=element_blank())+
      scale_colour_discrete(labels=c("Ajuste log-normal"))
    
    #Gráfico de la desviación estándar: (el gráfico coincide! pero el resultado NO!!)
    
    # plot(x=SDbinmedio, y=frecuencia, col = "red", main = "Distribución Log-Normal de SD", xlab = "Desviación estándar (m/s)", ylab = "Frecuencia", pch=2, lwd=2, ylim = c(0,(max(frecuencia1)+0.05)))
    # legend("topright", col = c("blue", "red"), legend = c("Ajuste log-normal", "Datos medidos"), lwd = 2, bty = "n", border = "black", box.lwd = 2, box.col = "black", box.lty = "l")
    # lines(x=SDbinmedio, y=frecuencia1, col="blue", lwd=2)
    
    Coeficientes<-list(mu=mu, sigma=sigma, Chi_Cuadrado=Chi_Cuadrado, RMSE=RMSE, R_2=R, graph=p)
    
    return(Coeficientes)
    
  }
  
  
  else{
    
    stop("method es un argumento inválido. Prueba 'MLE' ó 'Mom' ")
    
  }
  
  
  
  
  
  
  
  
}
