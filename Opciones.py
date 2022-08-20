"""
Created on Tuesday May 24 14:10:52 2022

@author: Osvaldo

opciones financieras
"""
import math
import copy
import tkinter as tk
import networkx as nx
from tkinter import *
from tkinter import ttk
from math import exp, factorial, log, sqrt
from matplotlib import pyplot as plt
import numpy as np 
import pandas as pd
import pandas_datareader.data as wb
from scipy.stats import norm
from PIL import Image,ImageTk


def MBS(q1,k,r,v,T):
    if v==0:
        data = pd.DataFrame()
        data = wb.DataReader(q1,'yahoo','2018-1-1')['Adj Close']
        log_returns = np.log(1+data.pct_change())
        S = data.iloc[-1]
        Vol = log_returns.std()*252**0.5
    else:
        S=q1
        Vol=v
    def d1(S,k,r,Vol,T):
        return (log(S/k)+((r + Vol**2)/2)*T)/(Vol * sqrt(T))

    def d2(S,k,r,Vol,T):
        return (log(S/k)+((r - Vol**2)/2)*T)/(Vol * sqrt(T))

    def BScall(S,k,r,Vol,T):
        d_uno = d1(S,k,r,Vol,T)
        d_dos = d2(S,k,r,Vol,T)
        return (S*norm.cdf(d_uno)) - (k*exp(-r*T)*norm.cdf(d_dos))

    def BSput(S,k,r,Vol,T):
        d_uno = d1(S,k,r,Vol,T)
        d_dos = d2(S,k,r,Vol,T)
        return (k*exp(-r*T)*norm.cdf(-d_dos)) - (S*norm.cdf(-d_uno)) 

    if v==0:
        print(data)
        print(Vol)

    tabla_viwBS(BScall(S,k,r,Vol,T),BSput(S,k,r,Vol,T),S,r,k,T,Vol)

def tabla_viwBS(cal1,put1,S0,r,k,T,v):
    
    main_window = tk.Tk()
    main_window.title("Tabla de valores")
    treeview = ttk.Treeview(main_window,columns=("size"))
    treeview.heading("#0", text="Descripción")
    treeview.heading("size", text="Valor")
    treeview.insert("",END,text="El precio de la acción es : ",values=(S0))
    treeview.insert("",END,text="Tasa de interes : ",values=(r))
    treeview.insert("",END,text="Precio de ejercicio: ",values=(k))
    treeview.insert("",END,text="Tiempo total (meses): ",values=(T))
    treeview.insert("",END,text="volatilidad : ",values=(v))
    treeview.insert("",END,text="El valor del Call es : ",values=(cal1))
    treeview.insert("",END,text="El valor del Put es : ",values=(put1))

    treeview.pack()
    main_window.mainloop()



def tabla_viw(S0,n,up,down,r,meses,k,call_calculado):
	
	main_window = tk.Tk()
	main_window.title("Tabla de valores")
	treeview = ttk.Treeview(main_window,columns=("size"))
	treeview.heading("#0", text="Archivo")
	treeview.heading("size", text="Tamaño")
	
	treeview.insert("",END,text="El precio de la acción es : ",values=(S0))
	treeview.insert("",END,text="Número de periodos : ",values=(n))
	treeview.insert("",END,text="Porcentaje que subio la acción : ",values=(up))
	treeview.insert("",END,text="Porcentaje que bajo la acción : ",values=(down))
	treeview.insert("",END,text="Tasa de interes : ",values=(r))
	treeview.insert("",END,text="Tiempo total : ",values=(meses))
	treeview.insert("",END,text="Precio de ejercicio : ",values=(k))
	treeview.insert("",END,text="El valor del Call es: ",values=(call_calculado))

	treeview.pack()
	main_window.mainloop()
	
def tasa(p,r,t,tt):
    if tt == 'c':
        F = p*exp(r*(t/12))
    else:
        F = p*(1 + (r/(12/tt)))**((t/tt))
    return round(F,8)

def binomial_vol(S0,n,T,vol,r,tt,K,opc,tipo):
    
    meses=T/n
    dt=T/(12*n)
    u=exp(sqrt(dt)*vol)
    d=1/u
    if r<log(u)/dt:
        print('Ta bien')
    else:
        print("No va a jalar")
    lc= [[S0]]
    ltemp =[]
    for i in range(n):
        for j,a in enumerate(lc[i]):
            if j == (len(lc[i])-1):
                ltemp.append(a*(u))
                ltemp.append(a*(d))
                lc.append(ltemp)
                ltemp=[]
            else:
                ltemp.append(a*(u))
    temp = tasa(1,-n*r*dt,12,tt)
    p=(exp(r*dt)-d)/(u-d)
    print('"prob"',round(p,8))

    lcall=[]
    lput=[]
    if opc=='c':
        lcall=copy.deepcopy(lc)
        for i in range(n,-1,-1):
            if i ==n:
                for j in range(len(lc[i])):
                    lcall[i][j] = round(max(lc[i][j]-K,0.0),8)
            elif i==0:
                for j in range(len(lc[i])):
                    lcall[i][j] = tasa(1,-r*dt,12,tt)*(lcall[i+1][j]*p + lcall[i+1][j+1]*(1-p))
            else:
                if tipo=="E":
                    for j in range(len(lc[i])):
                        lcall[i][j] = tasa(1,-r*dt,12,tt)*(lcall[i+1][j]*p + lcall[i+1][j+1]*(1-p))
                else:
                    for j in range(len(lc[i])):
                        lcall[i][j] = max(tasa(1,-r*dt,12,tt)*(lcall[i+1][j]*p + lcall[i+1][j+1]*(1-p)),lc[i][j]-K)
                        if n<20:
                            if lcall[i][j]==lc[i][j]-K:
                                print(f"Es optimo ejercitar el call en el nodo ubicado en {i},{len(lc[i])-1-2*j}") ##El nodo es optimo
        print('call',lcall[0][0])
    elif opc=='p':
        lput=copy.deepcopy(lc)
        for i in range(n,-1,-1):
            if i ==n:
                for j in range(len(lc[i])):
                    lput[i][j] = round(max(K-lc[i][j],0.0),8)
            elif i==0:
                for j in range(len(lc[i])):
                    lput[i][j] = tasa(1,-r*dt,12,tt)*(lput[i+1][j]*p + lput[i+1][j+1]*(1-p))
            else:
                if tipo=="E":
                    for j in range(len(lc[i])):
                        lput[i][j] = tasa(1,-r*dt,12,tt)*(lput[i+1][j]*p + lput[i+1][j+1]*(1-p))
                else:
                    for j in range(len(lc[i])):
                        lput[i][j] = max(tasa(1,-r*dt,12,tt)*(lput[i+1][j]*p + lput[i+1][j+1]*(1-p)),K-lc[i][j])
                        if n<20:
                            if lput[i][j]==K-lc[i][j]:
                                print(f"Es optimo ejercitar el put en el nodo ubicado en {i},{len(lc[i])-1-2*j}")##El nodo el optimo
        print('put',lput[0][0])
    else:
        lcall=copy.deepcopy(lc)
        for i in range(n,-1,-1):
            if i ==n:
                for j in range(len(lc[i])):
                    lcall[i][j] = round(max(lc[i][j]-K,0.0),8)
            elif i==0:
                for j in range(len(lc[i])):
                    lcall[i][j] = tasa(1,-r*dt,12,tt)*(lcall[i+1][j]*p + lcall[i+1][j+1]*(1-p))
            else:
                if tipo=="E":
                    for j in range(len(lc[i])):
                        lcall[i][j] = tasa(1,-r*dt,12,tt)*(lcall[i+1][j]*p + lcall[i+1][j+1]*(1-p))
                else:
                    for j in range(len(lc[i])):
                        lcall[i][j] = max(tasa(1,-r*dt,12,tt)*(lcall[i+1][j]*p + lcall[i+1][j+1]*(1-p)),lc[i][j]-K)
                        if n<20:
                            if lcall[i][j]==lc[i][j]-K:
                                print(f"Es optimo ejercitar el call en el nodo ubicado en {i},{len(lc[i])-1-2*j}") ##El nodo es optimo
        print('call',lcall[0][0])
        
        lput=copy.deepcopy(lc)
        for i in range(n,-1,-1):
            if i ==n:
                for j in range(len(lc[i])):
                    lput[i][j] = round(max(K-lc[i][j],0.0),8)
            elif i==0:
                for j in range(len(lc[i])):
                    lput[i][j] = tasa(1,-r*dt,12,tt)*(lput[i+1][j]*p + lput[i+1][j+1]*(1-p))
            else:
                if tipo=="E":
                    for j in range(len(lc[i])):
                        lput[i][j] = tasa(1,-r*dt,12,tt)*(lput[i+1][j]*p + lput[i+1][j+1]*(1-p))
                else:
                    for j in range(len(lc[i])):
                        lput[i][j] = max(tasa(1,-r*dt,12,tt)*(lput[i+1][j]*p + lput[i+1][j+1]*(1-p)),K-lc[i][j])
                        if n<20:
                            if lput[i][j]==K-lc[i][j]:
                                print(f"Es optimo ejercitar el put en el nodo ubicado en {i},{len(lc[i])-1-2*j}")##El nodo el optimo
        print('put',lput[0][0])
        Arbitraje_forwd(lput[0][0], lcall[0][0], K, T, S0, r, tt)
    if n<50:
        grafos(lc,lcall,lput,opc,n)
    #grafos(lc,lcall,lput,opc,n)
    
def grafos(lc,lcall,lput,opc,n):
    if opc=='c':
        G=nx.Graph()
        l=copy.deepcopy(lc)
        mm=0
        for i in range(len(lc)):
            for j in range(len(lc[i])):
                l[i][j]=mm
                mm+=1
        for i in range(len(lc)):
            nn=len(lc[i])
            k=0
            for j in range(nn-1,-nn,-(2)):
                G.add_node(l[i][k],precio = round(lc[i][k],4),call = round(lcall[i][k],4),pos=(i,j))
                if i>0 and k<i:
                    G.add_edge(l[i-1][k],l[i][k])
                    G.add_edge(l[i-1][k],l[i][k+1])
                k+=1
        pos = nx.get_node_attributes(G,'pos')
        nl=nx.get_node_attributes(G,'precio')
        nl1=nx.get_node_attributes(G,'call')
        print("arbol de precios")
        if n<=10:
            print(lc)
        plt.figure(figsize=(100, 100))
        plt.axis('off')
        nx.draw(G,pos,node_size=2000/n)
        nx.draw_networkx_labels(G,pos,labels=nl)
        plt.show()
        print("arbol de precios del call")
        if n<=10:
            print(lcall)
        plt.figure(figsize=(100, 100))
        plt.axis('off')
        nx.draw(G,pos,node_size=2000/n)
        nx.draw_networkx_labels(G,pos,labels=nl1)
        plt.show()
    elif opc=='p':
        G=nx.Graph()
        l=copy.deepcopy(lc)
        mm=0
        for i in range(len(lc)):
            for j in range(len(lc[i])):
                l[i][j]=mm
                mm+=1
        for i in range(len(lc)):
            nn=len(lc[i])
            k=0
            for j in range(nn-1,-nn,-(2)):
                G.add_node(l[i][k],precio = round(lc[i][k],4),put = round(lput[i][k],4),pos=(i,j))
                if i>0 and k<i:
                    G.add_edge(l[i-1][k],l[i][k])
                    G.add_edge(l[i-1][k],l[i][k+1])
                k+=1
        pos = nx.get_node_attributes(G,'pos')
        nl=nx.get_node_attributes(G,'precio')
        nl2=nx.get_node_attributes(G,'put')
        print("arbol de precios")
        if n<=10:
            print(lc)
        plt.figure(figsize=(100, 100))
        plt.axis('off')
        nx.draw(G,pos,node_size=2000/n)
        nx.draw_networkx_labels(G,pos,labels=nl)
        plt.show()
        print("arbol de precios del put")
        if n<=10:
            print(lput)
        plt.figure(figsize=(100, 100))
        plt.axis('off')
        nx.draw(G,pos,node_size=2000/n)
        nx.draw_networkx_labels(G,pos,labels=nl2)
        plt.show()
    else:
        G=nx.Graph()
        l=copy.deepcopy(lc)
        mm=0
        for i in range(len(lc)):
            for j in range(len(lc[i])):
                l[i][j]=mm
                mm+=1
        for i in range(len(lc)):
            nn=len(lc[i])
            k=0
            for j in range(nn-1,-nn,-(2)):
                G.add_node(l[i][k],precio = round(lc[i][k],4),call = round(lcall[i][k],4),put = round(lput[i][k],4),pos=(i,j))
                if i>0 and k<i:
                    G.add_edge(l[i-1][k],l[i][k])
                    G.add_edge(l[i-1][k],l[i][k+1])
                k+=1
        pos = nx.get_node_attributes(G,'pos')
        nl=nx.get_node_attributes(G,'precio')
        nl1=nx.get_node_attributes(G,'call')
        nl2=nx.get_node_attributes(G,'put')
        print("arbol de precios")
        if n<=10:
            print(lc)
        plt.figure(figsize=(100,100))
        plt.axis('off')
        nx.draw(G,pos,node_size=2000/n)
        nx.draw_networkx_labels(G,pos,labels=nl)
        plt.show()
        print("arbol de precios del call")
        if n<=10:
            print(lcall)
        plt.figure(figsize=(100, 100))
        plt.axis('off')
        nx.draw(G,pos,node_size=2000/n)
        nx.draw_networkx_labels(G,pos,labels=nl1)
        plt.show()
        print("arbol de precios del put")
        if n<=10:
            print(lput)
        plt.figure(figsize=(100, 100))
        plt.axis('off')
        nx.draw(G,pos,node_size=2000/n)
        nx.draw_networkx_labels(G,pos,labels=nl2)
        plt.show()

def Arbitraje_forwd(p, c, K, t, S0, r, tt):
    
    print('\n Recordemos que la condición para o tener arbitraje es:')
    print('c+Ke^{-rt} = p + S0')
    print(f'Al despejar "p" y "c" de la ecuacion respectivamente obtenemos:')
    p1=c+tasa(K,-r,t,tt) - S0
    print(f'El valor del put deberia ser {p1}')
    c1 = p +S0 -tasa(K,-r,t,tt)#K*exp(-r*t/12)
    print(f'El valor del call deberia ser {c1}')
    if round(p,4)==round(p1,4) and round(c,4)==round(c1,4):
        print('Por lo que no existe oportunidad de arbitraje \n')
    else:
        print('Por lo que si existe oportunidad de arbitraje \n')

class SampleApp(tk.Tk):
    def __init__(self):
        Tk.__init__(self)
        self._frame = None
        self.switch_frame(StartPage)

    def switch_frame(self, frame_class):
        new_frame = frame_class(self)
        if self._frame is not None:
            self._frame.destroy()
        self._frame = new_frame
        new_frame.config(width="400",height="480")
        new_frame.pack(padx=5,pady=5,ipadx=5,ipady=5)
        self._frame.pack()

class StartPage(tk.Frame):
    def __init__(self, master):
        Frame.__init__(self, master)
        Label(self, text="Seleccione el modelo que desea utilizar").pack(side="top", fill="x", pady=10)
        Button(self, text="Binomial",command=lambda: master.switch_frame(ModeloBinomial)).pack(padx=5,pady=5,ipadx=2.5,ipady=2.5)
        Button(self, text="Black - Scholes",command=lambda: master.switch_frame(SeleccionVolatilidadBS)).pack(padx=5,pady=5,ipadx=2.5,ipady=2.5)

class SeleccionVolatilidadBS(tk.Frame):
    def __init__(self,master):
        Frame.__init__(self,master)
        Label(self,text="Desea ingresar la volatilidad").pack(padx=5,pady=5,ipadx=2.5,ipady=2.5)
        Button(self,text="Click here!",command=lambda: master.switch_frame(ModeloBlackScholesCV)).pack(padx=5,pady=5,ipadx=2.5,ipady=2.5)
        Label(self,text="Desea calcular la volatilidad de una acción").pack(padx=5,pady=5,ipadx=2.5,ipady=2.5)
        Button(self,text="Click here!",command=lambda: master.switch_frame(ModeloBlackScholes)).pack(padx=5,pady=5,ipadx=2.5,ipady=2.5)

class ModeloBlackScholesCV(tk.Frame):
    def __init__(self,master):
        Frame.__init__(self,master)
        def datos_guardados():
            Ac_data = float(str(Ac.get()))
            ts_data = float(str(ts.get()))*0.01
            ejer_data = int(str(ejer.get()))
            tim_data = float(str(tim.get()))/12
            vol_data = float(str(vol.get()))*0.01

            MBS(Ac_data,ejer_data,ts_data,vol_data,tim_data)

        Label(self,text="Black - Scholes").pack()
        app.geometry('400x480')
        Label(text="Ingrese el precio de la acción actual").place(x=22,y=20) 
        Label(text="Ingrese la tasa libre de riesgo en porcentaje (r)").place(x=22,y=70)
        Label(text="Ingrese el precio de ejercicio (k)").place(x=22,y=110)
        Label(text="Ingrese el tiempo en meses (T)").place(x=22,y=170)
        Label(text="Ingrese la volatilidad en porcentaje").place(x=22,y=230)
        Button(app,text="Enter",command=datos_guardados).place(x=22,y=299)
        Button(app,text="Home",command=lambda: master.switch_frame(StartPage)).place(x=92,y=299)

        Ac = StringVar()
        ts = StringVar()
        ejer = StringVar()
        tim = StringVar()
        vol = StringVar()

        Ac_entry = Entry(textvariable=Ac,width="40")
        ts_entry = Entry(textvariable=ts,width="40")
        ejer_entry = Entry(textvariable=ejer,width="40")
        tim_entry = Entry(textvariable=tim,width="40")
        vol_entry  = Entry(textvariable=vol,width="40")

        Ac_entry.place(x=22,y=40)
        ts_entry.place(x=22,y=90)
        ejer_entry.place(x=22,y=150)
        vol_entry.place(x=22,y=265)
        tim_entry.place(x=22,y=199)

class ModeloBlackScholes(tk.Frame):
    def __init__(self,master):
        Frame.__init__(self,master)
        def datos_guardados():
            Ac_data = str(Ac.get())
            ts_data = float(str(ts.get()))*0.01
            ejer_data = int(str(ejer.get()))
            tim_data = float(str(tim.get()))/12

            MBS(Ac_data,ejer_data,ts_data,0,tim_data)
    
        Label(self,text="Black - Scholes").pack()
        app.geometry('400x480')
        Label(text="Ingrese la acción que desea cotizar").place(x=22,y=20)
        Label(text="Ingrese la tasa libre de riesgo en porcentaje (r)").place(x=22,y=70)
        Label(text="Ingrese el precio de ejercicio (k)").place(x=22,y=110)
        Label(text="Ingrese el tiempo en meses (T)").place(x=22,y=170)
        Button(app,text="Enter",command=datos_guardados).place(x=22,y=270)
        Button(app,text="Home", 
          command=lambda: master.switch_frame(StartPage)).place(x=92,y=270)
        Ac = StringVar()
        ts = StringVar()
        ejer = StringVar()
        tim = StringVar()

        Ac_entry = Entry(textvariable=Ac,width="40")
        ts_entry = Entry(textvariable=ts,width="40")
        ejer_entry = Entry(textvariable=ejer,width="40")
        tim_entry = Entry(textvariable=tim,width="40")

        Ac_entry.place(x=22,y=40)
        ts_entry.place(x=22,y=90)
        ejer_entry.place(x=22,y=150)
        tim_entry.place(x=22,y=199)

class ModeloBinomial(tk.Frame):
    def __init__(self,master):
        Frame.__init__(self,master)
        Label(self, text="Seleccione su tipo de opción financiera").pack(side="top",fill="x",pady=10)
        Button(self, text="Americana",
                  command=lambda: master.switch_frame(Ingresa_datosAm)).pack(padx=5,pady=5,ipadx=2.5,ipady=2.5)
        Button(self, text="Europera",
                  command=lambda: master.switch_frame(Ingresa_datosEu)).pack(padx=5,pady=5,ipadx=2.5,ipady=2.5)

class Ingresa_datosAm(tk.Frame):
	def __init__(self,master):
		Frame.__init__(self,master)
		def datos_guardadosCall():

			So_data = int(str(So.get()))
			n_data = int(str(n.get()))
			v_data = float(str(v.get()))*0.01
			r_data = float(str(r.get()))*0.01
			T_data = float(str(T.get()))
			k_data = int(str(k.get()))
			
			binomial_vol(So_data,n_data,T_data,v_data,r_data,'c',k_data,'c','A')
		def datos_guardadosPut():

			So_data = int(str(So.get()))
			n_data = int(str(n.get()))
			v_data = float(str(v.get()))*0.01
			r_data = float(str(r.get()))*0.01
			T_data = float(str(T.get()))
			k_data = int(str(k.get()))
			
			binomial_vol(So_data,n_data,T_data,v_data,r_data,'c',k_data,'p','A')
		
		app.geometry('400x480')

		Label(text="Ingrese el precio de la acción actual (So)").place(x=22,y=20)
		Label(text="Ingrese el número de periodos (n)").place(x=22,y=70)
		Label(text="Ingrese la volatilidad de So en porcentaje (v)").place(x=22,y=140)
		Label(text="Ingrese la tasa de interes en porcentaje (r)").place(x=22,y=200)
		Label(text="Ingrese el tiempo total en meses (T)").place(x=22,y=260)
		Label(text="Ingrese el precio de ejercicio (k)").place(x=22,y=310)
		Label(text="¿Que desea calcular?").place(x=22,y=370)
		Button(app,text="Call",command=datos_guardadosCall).place(x=22,y=400)
		Button(app,text="put",command=datos_guardadosPut).place(x=62,y=400)
		Button(self, text="Home",
                  command=lambda: master.switch_frame(StartPage)).place(x=302,y=130)

		So = StringVar()
		n = StringVar()
		v = StringVar()
		r = StringVar()
		T = StringVar()
		k = StringVar()

		So_entry = Entry(textvariable=So,width="40")
		n_entry = Entry(textvariable=n,width="40")
		v_entry = Entry(textvariable=v,width="40")
		r_entry = Entry(textvariable=r,width="40")
		T_entry = Entry(textvariable=T,width="40")
		k_entry = Entry(textvariable=k,width="40")

		So_entry.place(x=22,y=40)
		n_entry.place(x=22,y=100)
		v_entry.place(x=22,y=170)
		r_entry.place(x=22,y=230)
		T_entry.place(x=22,y=290)
		k_entry.place(x=22,y=340)

class Ingresa_datosEu(tk.Frame):
	def __init__(self,master):
		Frame.__init__(self,master)
		def datos_guardadosCall():

			So_data = int(str(So.get()))
			n_data = int(str(n.get()))
			v_data = float(str(v.get()))*0.01
			r_data = float(str(r.get()))*0.01
			T_data = float(str(T.get()))
			k_data = int(str(k.get()))
			
			binomial_vol(So_data,n_data,T_data,v_data,r_data,'c',k_data,'c','E')
		def datos_guardadosPut():

			So_data = int(str(So.get()))
			n_data = int(str(n.get()))
			v_data = float(str(v.get()))*0.01
			r_data = float(str(r.get()))*0.01
			T_data = float(str(T.get()))
			k_data = int(str(k.get()))
			
			binomial_vol(So_data,n_data,T_data,v_data,r_data,'c',k_data,'p','E')
		
		app.geometry('400x480')

		Label(text="Ingrese el precio de la acción actual (So)").place(x=22,y=20)
		Label(text="Ingrese el número de periodos (n)").place(x=22,y=70)
		Label(text="Ingrese la volatilidad de So en porcentaje (v)").place(x=22,y=140)
		Label(text="Ingrese la tasa de interes en porcentaje (r)").place(x=22,y=200)
		Label(text="Ingrese el tiempo total en meses (T)").place(x=22,y=260)
		Label(text="Ingrese el precio de ejercicio (k)").place(x=22,y=310)
		Label(text="¿Que desea calcular?").place(x=22,y=370)
		Button(app,text="Call",command=datos_guardadosCall).place(x=22,y=400)
		Button(app,text="put",command=datos_guardadosPut).place(x=62,y=400)
		Button(self, text="Home",
                  command=lambda: master.switch_frame(StartPage)).place(x=302,y=130)

		So = StringVar()
		n = StringVar()
		v = StringVar()
		r = StringVar()
		T = StringVar()
		k = StringVar()

		So_entry = Entry(textvariable=So,width="40")
		n_entry = Entry(textvariable=n,width="40")
		v_entry = Entry(textvariable=v,width="40")
		r_entry = Entry(textvariable=r,width="40")
		T_entry = Entry(textvariable=T,width="40")
		k_entry = Entry(textvariable=k,width="40")

		So_entry.place(x=22,y=40)
		n_entry.place(x=22,y=100)
		v_entry.place(x=22,y=170)
		r_entry.place(x=22,y=230)
		T_entry.place(x=22,y=290)
		k_entry.place(x=22,y=340)

if __name__ == "__main__":
    app = SampleApp()
    app.geometry('400x480')
    app.mainloop()
