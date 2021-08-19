from RungeKutta4th import RungeKutta4th
import numpy as np
import tkinter as tk
import time
from tkinter import font
from math import *


def set_m1(event):
    button_start["state"] = "disabled"
    label_m1_value["text"] = f"{10**(var_m1.get()):.2f}"
    update_pendulum(var_m1.get(), var_m2.get(), 10**var_l1.get(), 10**var_l2.get(), var_t1.get(), var_t2.get())

def set_m2(event):
    button_start["state"] = "disabled"
    label_m2_value["text"] = f"{10**(var_m2.get()):.2f}"
    update_pendulum(var_m1.get(), var_m2.get(), 10**var_l1.get(), 10**var_l2.get(), var_t1.get(), var_t2.get())

def set_l1(event):
    button_start["state"] = "disabled"
    label_l1_value["text"] = f"{10**(var_l1.get()):.2f}"
    update_pendulum(var_m1.get(), var_m2.get(), 10**var_l1.get(), 10**var_l2.get(), var_t1.get(), var_t2.get())

def set_l2(event):
    button_start["state"] = "disabled"
    label_l2_value["text"] = f"{10**(var_l2.get()):.2f}"
    update_pendulum(var_m1.get(), var_m2.get(), 10**var_l1.get(), 10**var_l2.get(), var_t1.get(), var_t2.get())

def set_g(event):
    button_start["state"] = "disabled"
    label_g_value["text"] = f"{9.80665*10**(var_g.get()):.2f}"

def set_t1(event):
    button_start["state"] = "disabled"
    label_t1_value["text"] = f"{var_t1.get()*180/pi:.1f}°"
    update_pendulum(var_m1.get(), var_m2.get(), 10**var_l1.get(), 10**var_l2.get(), var_t1.get(), var_t2.get())

def set_t2(event):
    button_start["state"] = "disabled"
    label_t2_value["text"] = f"{var_t2.get()*180/pi:.1f}°"
    update_pendulum(var_m1.get(), var_m2.get(), 10**var_l1.get(), 10**var_l2.get(), var_t1.get(), var_t2.get())

def set_default():
    button_start["state"] = "disabled"
    var_m1.set(0)
    var_m2.set(0)
    var_l1.set(0)
    var_l2.set(0)
    var_g.set(0)
    var_t1.set(0)
    var_t2.set(0)
    label_m1_value["text"] = f"{10**(var_m1.get()):.2f}"
    label_m2_value["text"] = f"{10**(var_m2.get()):.2f}"
    label_l1_value["text"] = f"{10**(var_l1.get()):.2f}"
    label_l2_value["text"] = f"{10**(var_l2.get()):.2f}"
    label_g_value["text"] = f"{9.80665*10**(var_g.get()):.2f}"
    label_t1_value["text"] = f"{var_t1.get()*180/pi:.1f}°"
    label_t2_value["text"] = f"{var_t2.get()*180/pi:.1f}°"
    update_pendulum(var_m1.get(), var_m2.get(), 10**var_l1.get(), 10**var_l2.get(), var_t1.get(), var_t2.get())

def update_pendulum(m1, m2, l1, l2, t1, t2):
    r_ball = 10
    r_stick = 175

    ball_r1, ball_r2 = int(r_ball*(1.1+m1)), int(r_ball*(1.1+m2))
    stick_r1, stick_r2 = int(l1*r_stick/(l1+l2)), int(l2*r_stick/(l1+l2))

    pos1_x, pos1_y = 200, 200
    pos2_x, pos2_y = pos1_x + stick_r1*sin(t1), pos1_y + stick_r1*cos(t1)
    pos3_x, pos3_y = pos2_x + stick_r2*sin(t2), pos2_y + stick_r2*cos(t2)

    canvas_pendulum.coords(line1, pos1_x, pos1_y, pos2_x, pos2_y)
    canvas_pendulum.coords(line2, pos2_x, pos2_y, pos3_x, pos3_y)
    canvas_pendulum.coords(ball1, pos2_x-ball_r1, pos2_y-ball_r1, pos2_x+ball_r1, pos2_y+ball_r1)
    canvas_pendulum.coords(ball2, pos3_x-ball_r2, pos3_y-ball_r2, pos3_x+ball_r2, pos3_y+ball_r2)

    label_pendulum_t1_val["text"] = f"{(t1*180/pi)%360:.1f}°"
    label_pendulum_t2_val["text"] = f"{(t2*180/pi)%360:.1f}°"

    window.update()



def update_consistently():
    global t1, t2, w1, w2, t
    m1, m2 = 10**var_m1.get(), 10**var_m2.get()
    l1, l2 = 10**var_l1.get(), 10**var_l2.get()
    g = 9.80665*10**var_g.get()

    n = 4
    f1 = lambda t, x: x[2]
    f2 = lambda t, x: x[3]
    def f3(t, x):
        numer1 = -g*(2*m1+m2)*np.sin(x[0]) -m2*g*np.sin(x[0]-2*x[1])
        numer2 = -2*np.sin(x[0]-x[1])*m2*(x[3]**2*l2+x[2]**2*l1*np.cos(x[0]-x[1]))
        numer = numer1 + numer2
        denom = l1*(2*m1+m2-m2*np.cos(2*x[0]-2*x[1]))
        return numer/denom
    def f4(t, x):
        numer1 = x[2]**2*l1*(m1+m2) + g*(m1+m2)*np.cos(x[0])
        numer2 = x[3]**2*l2*m2*np.cos(x[0]-x[1])
        numer = 2*np.sin(x[0]-x[1])*(numer1+numer2)
        denom = l2*(2*m1+m2-m2*np.cos(2*x[0]-2*x[1]))
        return numer/denom
    f = np.array([f1, f2, f3, f4])
    x = np.array([t1, t2, w1, w2])

    fps = 50
    step = 10
    dt = 1/fps

    while state == 1:
        start_time = time.time()
        
        t, x = RungeKutta4th(n, f, t, x, t+dt, step)
        t1, t2, w1, w2 = x

        update_pendulum(log10(m1), log10(m2), l1, l2, t1, t2)
        label_pendulum_t_val["text"] = f"{t:.1f}s"

        end_time = time.time()
        elapsed = end_time-start_time
        time.sleep(max(dt-elapsed, 0))



def get_ready():
    global t
    global state
    t = 0
    state = 0
    button_set["text"] = "set"
    button_start["text"] = "start"
    button_start["state"] = "normal"
    button_pause["state"] = "disabled"
    button_set_default["state"] = "normal"
    scale_m1["state"] = "normal"
    scale_m2["state"] = "normal"
    scale_l1["state"] = "normal"
    scale_l2["state"] = "normal"
    scale_g["state"] = "normal"
    scale_t1["state"] = "normal"
    scale_t2["state"] = "normal"
    update_pendulum(var_m1.get(), var_m2.get(), 10**var_l1.get(), 10**var_l2.get(), var_t1.get(), var_t2.get())



def start():
    global t1, t2, w1, w2, t
    global state
    t1, t2 = var_t1.get(), var_t2.get()
    w1, w2 = 0, 0
    t = 0
    state = 1

    button_set["text"] = "reset"
    button_start["text"] = "restart"
    button_pause["state"] = "normal"
    button_continue["state"] = "disabled"
    button_set_default["state"] = "disabled"
    scale_m1["state"] = "disabled"
    scale_m2["state"] = "disabled"
    scale_l1["state"] = "disabled"
    scale_l2["state"] = "disabled"
    scale_g["state"] = "disabled"
    scale_t1["state"] = "disabled"
    scale_t2["state"] = "disabled"

    update_consistently()


def pause():
    global state
    state = 2
    button_pause["state"] = "disabled"
    button_continue["state"] = "normal"


def continue_():
    global t1, t2
    global state
    state = 1
    button_pause["state"] = "normal"
    button_continue["state"] = "disabled"

    update_consistently()





state = 0 # 0: before start // 1: ongoing // 2: paused
t1, t2 = 0, 0
w1, w2 = 0, 0
t = 0

window = tk.Tk()
window.title("Double Pendulum Modeling")
window.geometry("820x570+0+0")
window.resizable(False, False)

var_m1 = tk.DoubleVar()
var_m2 = tk.DoubleVar()
var_l1 = tk.DoubleVar()
var_l2 = tk.DoubleVar()
var_g = tk.DoubleVar()
var_t1 = tk.DoubleVar()
var_t2 = tk.DoubleVar()

font_title = font.Font(family="Arial", size=20, weight="bold")
font_int_big = font.Font(family="Arial", size=15)
font_int_small = font.Font(family="Arial", size=12)

label_title = tk.Label(window, text="Double Pendulum", font=font_title)
label_title.place(x=10, y=10, width=800, height=50)

canvas_pendulum = tk.Canvas(window, relief="solid", bd=1, bg="white")
canvas_pendulum.place(x=10, y=60, width=400, height=400)
line1 = canvas_pendulum.create_line(200, 200, 200, 287, width=3)
line2 = canvas_pendulum.create_line(200, 287, 200, 375, width=3)
ball1 = canvas_pendulum.create_oval(189, 276, 211, 298, fill="blue", width=3)
ball2 = canvas_pendulum.create_oval(189, 364, 211, 386, fill="blue", width=3)

#canvas_graph = tk.Canvas(window, relief="solid", bd=1, bg="white")
#canvas_graph.place(x=10, y=460, width=400, height=250)
label_pendulum_t = tk.Label(window, text="t=", font=font_int_small)
label_pendulum_t.place(x=10, y=460, width=30, height=50)
label_pendulum_t_val = tk.Label(window, text="0.0s", font=font_int_small, anchor="e")
label_pendulum_t_val.place(x=40, y=460, width=50, height=50)

label_pendulum_t1 = tk.Label(window, text="θ1=", font=font_int_small)
label_pendulum_t1.place(x=130, y=460, width=30, height=50)
label_pendulum_t1_val = tk.Label(window, text="0.0°", font=font_int_small, anchor="e")
label_pendulum_t1_val.place(x=160, y=460, width=50, height=50)

label_pendulum_t2 = tk.Label(window, text="θ2=", font=font_int_small)
label_pendulum_t2.place(x=250, y=460, width=30, height=50)
label_pendulum_t2_val = tk.Label(window, text="0.0°", font=font_int_small, anchor="e")
label_pendulum_t2_val.place(x=280, y=460, width=50, height=50)


label_m1 = tk.Label(window, text="m1", font=font_int_big)
label_m1.place(x=410, y=100, width=100, height=50)
scale_m1 = tk.Scale(window, from_=-1, to=1, resolution=0.01, orient="horizontal", showvalue=False, command=set_m1, variable=var_m1)
scale_m1.place(x=510, y=120, width=200, height=50)
label_m1_value = tk.Label(window, text="1.00", font=font_int_big)
label_m1_value.place(x=710, y=100, width=100, height=50)

label_m2 = tk.Label(window, text="m2", font=font_int_big)
label_m2.place(x=410, y=150, width=100, height=50)
scale_m2 = tk.Scale(window, from_=-1, to=1, resolution=0.01, orient="horizontal", showvalue=False, command=set_m2, variable=var_m2)
scale_m2.place(x=510, y=170, width=200, height=50)
label_m2_value = tk.Label(window, text="1.00", font=font_int_big)
label_m2_value.place(x=710, y=150, width=100, height=50)

label_l1 = tk.Label(window, text="l1", font=font_int_big)
label_l1.place(x=410, y=200, width=100, height=50)
scale_l1 = tk.Scale(window, from_=-1, to=1, resolution=0.01, orient="horizontal", showvalue=False, command=set_l1, variable=var_l1)
scale_l1.place(x=510, y=220, width=200, height=50)
label_l1_value = tk.Label(window, text="1.00", font=font_int_big)
label_l1_value.place(x=710, y=200, width=100, height=50)

label_l2 = tk.Label(window, text="l2", font=font_int_big)
label_l2.place(x=410, y=250, width=100, height=50)
scale_l2 = tk.Scale(window, from_=-1, to=1, resolution=0.01, orient="horizontal", showvalue=False, command=set_l2, variable=var_l2)
scale_l2.place(x=510, y=270, width=200, height=50)
label_l2_value = tk.Label(window, text="1.00", font=font_int_big)
label_l2_value.place(x=710, y=250, width=100, height=50)

label_g = tk.Label(window, text="g", font=font_int_big)
label_g.place(x=410, y=300, width=100, height=50)
scale_g = tk.Scale(window, from_=-1, to=1, resolution=0.01, orient="horizontal", showvalue=False, command=set_g, variable=var_g)
scale_g.place(x=510, y=320, width=200, height=50)
label_g_value = tk.Label(window, text="9.81", font=font_int_big)
label_g_value.place(x=710, y=300, width=100, height=50)

label_t1 = tk.Label(window, text="θ1", font=font_int_big)
label_t1.place(x=410, y=350, width=100, height=50)
scale_t1 = tk.Scale(window, from_=0, to=2*pi, resolution=0.01, orient="horizontal", showvalue=False, command=set_t1, variable=var_t1)
scale_t1.place(x=510, y=370, width=200, height=50)
label_t1_value = tk.Label(window, text="0.0°", font=font_int_big)
label_t1_value.place(x=710, y=350, width=100, height=50)

label_t2 = tk.Label(window, text="θ2", font=font_int_big)
label_t2.place(x=410, y=400, width=100, height=50)
scale_t2 = tk.Scale(window, from_=0, to=2*pi, resolution=0.01, orient="horizontal", showvalue=False, command=set_t2, variable=var_t2)
scale_t2.place(x=510, y=420, width=200, height=50)
label_t2_value = tk.Label(window, text="0.0°", font=font_int_big)
label_t2_value.place(x=710, y=400, width=100, height=50)

button_set_default = tk.Button(window, text="set to default...", font=font_int_small, command=set_default)
button_set_default.place(x=610, y=460, width=200, height=30)

button_set = tk.Button(window, text="set", font=font_int_big, command=get_ready)
button_set.place(x=410, y=510, width=100, height=50)
button_start = tk.Button(window, text="start", font=font_int_big, command=start, state="disabled")
button_start.place(x=510, y=510, width=100, height=50)
button_pause = tk.Button(window, text="pause", font=font_int_big, command=pause, state="disabled")
button_pause.place(x=610, y=510, width=100, height=50)
button_continue = tk.Button(window, text="continue", font=font_int_big, command=continue_, state="disabled")
button_continue.place(x=710, y=510, width=100, height=50)

set_default()

window.mainloop()
