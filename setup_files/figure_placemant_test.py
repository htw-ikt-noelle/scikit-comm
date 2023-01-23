import matplotlib.pyplot as plt
import tkinter as tk

plt.close('all')

for nf in range(11):
    plt.figure()

# 
taskbar_offset = 1 * 40
figure_toolbar = 1 * 64
# offset from top left screen edge [x, y]
offset = [1*1920, 0*100]
screen_resolution = [1920, 1080]

nc = 4
nr = 3

auto_layout = False

figHandle = list(map(plt.figure, plt.get_fignums()))   
n_fig = len(figHandle)

if n_fig <= 0:
    raise ValueError('no figures found to place')

if screen_resolution:
    screen_resolution[1] = screen_resolution[1] -  taskbar_offset
else:    
    # workaround to determine screen resolution:
    # * open tk window
    # * place the window according to offset
    # * make it fullscreen
    # * get width and height
    # * kill window
    # from https://stackoverflow.com/questions/3129322/how-do-i-get-monitor-resolution-in-python
    root = tk.Tk()
    root.update_idletasks()
    root.geometry(f'100x100+{offset[0]}+{offset[1]}')
    # print(root.geometry())
    root.attributes('-fullscreen', True)
    # print(root.geometry())
    root.state('iconic')
    screen_resolution = []
    screen_resolution.append(root.winfo_screenwidth())
    # reduce scrren height by the (windows) task bar
    screen_resolution.append(root.winfo_screenheight() -  taskbar_offset)
    root.destroy()

# auto layout?
if auto_layout:
    grid = [
        [1,1],[1,2],
        [2,2],[2,2],
        [2,3],[2,3],
        [3,3],[3,3],[3,3],
        [3,4],[3,4],[3,4],
        [4,4],[4,4],[4,4],[4,4],
        [4,5],[4,5],[4,5],[4,5],
        [4,6],[4,6],[4,6],[4,6],
        [4,7],[4,7],[4,7],[4,7],
        [4,8],[4,8],[4,8],[4,8]
        ]
   
    if n_fig > len(grid)*2:
        raise ValueError('more figures opened than layout options available')        
    
    # portrait mode
    if screen_resolution[0] < screen_resolution[1]:
        nc = grid[n_fig-1][0]
        nr = grid[n_fig-1][1]
    # landscape mode
    else:
        nc = grid[n_fig-1][1]
        nr = grid[n_fig-1][0]
# manual layout
else:
    if (nc * nr) < n_fig:
        raise ValueError(f'more figures opened ({n_fig}) than rows times coloumns given ({nc*nr}): try to increase numbers or switch to auto layout mode')
        

fig_width = screen_resolution[0]/nc 
fig_height = screen_resolution[1]/nr - figure_toolbar 

fig_cnt = 0
for r in range(nr):
    for c in range(nc):
        if fig_cnt >= n_fig:
            break        
        figHandle[fig_cnt].set_figheight(fig_height / figHandle[fig_cnt].get_dpi())
        figHandle[fig_cnt].set_figwidth(fig_width / figHandle[fig_cnt].get_dpi())
        if r == 0:
            figHandle[fig_cnt].canvas.manager.window.move(int(fig_width*c + offset[0]), int(fig_height*r) + offset[1])
        else:
            figHandle[fig_cnt].canvas.manager.window.move(int(fig_width*c + offset[0]), int((fig_height+figure_toolbar)*r) + offset[1])
        fig_cnt += 1