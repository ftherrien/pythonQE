from copy import deepcopy
import numpy as np

tol = 1e-10

def reciprocal(cell):
    rec = deepcopy(cell)
    V = np.linalg.det(cell)
    for i in range(3):
        rec[i,:] = 2*np.pi/V*np.cross(cell[:,(i+1)%3], cell[:,(i+2)%3])
    return rec

def explicit_path(path):
    epath = []
    for i,line in enumerate(path[:-1].T):
        for j in range(len(line))[:-1]:
            if j == 0:
                linpath = np.reshape(np.linspace(line[j],path[j,i+1],int(line[-1])+1), (int(line[-1])+1),1)
            else:
                linpath = np.concatenate((linpath, np.reshape(np.linspace(line[j],path[j,i+1],int(line[-1])+1), (int(line[-1])+1,1))), axis=1)
        linpath = np.concatenate((linpath, np.ones((int(line[-1])+1),1)),axis=1)
        if i == 0:
            epath = linpath 
        else:
            epath = np.concatenate((epath, linpath), axis = 0)
    return epath

def unique(closest):
    unique = []
    for line in closest.T:
        there = False
        for check in unique:
            if all(check == line):
                there = True
                break
        if not there:
            unique.append(line)
    return np.array(unique).T

def approx(i, vec):
    if i:
        return np.ceil(vec-tol)
    else:
        return np.floor(vec+tol)

def closest_box(path, rsc, irsc):
    epath = explicit_path(path)*np.pi*2 # explicit_path
    epathrc = irsc.dot(epath[0:3,:]) # explicit_reciprocal_path

    closest = [] 
    for i in range(2):
        for j in range(2):
            for k in range(2):
                closest.append(np.concatenate((approx(i,epathrc[0:1,:]), approx(j,epathrc[1:2,:]), approx(k,epathrc[2:3, :])), axis=0))

    newpath = unique(rsc.dot(np.concatenate(closest, axis=1)))/ (np.pi*2)

    return np.concatenate((newpath, np.ones((1,np.shape(newpath)[0]))), axis=0) #Adding ones for QE

def on_path(path, rsc, irsc):
    epath = explicit_path(path)*np.pi*2
    epathrc = irsc.dot(epath[0:3,:])
    
    npts = path[3,:]
    path = path[:3,:]

    closest_points = rsc.dot(unique(np.round(epathrc))) / (2*np.pi)
    onpath = []
    for i,v in enumerate(path[:,:-1].T):
        if npts[i] > 1:
            for p in closest_points.T:
                isonpath = True
                t = []
                for j in range(3):
                    if abs((path[j, i+1]-v[j])) >= tol: # make sure there is no division by zero
                        t.append((p[j]-v[j])/(path[j, i+1]-v[j]))
                    else:
                        if abs((p[j]-v[j])) >= tol: # if so, check if nominator is near 0
                            isonpath = False
            
                # Makes sure the multiplier t is the same for each component
                if len(t) == 2:
                    if abs(t[1] - t[0]) >= tol:
                        isonpath = False
                elif len(t) == 3:
                    for j in range(3):
                        if abs(t[j]-t[(j+1)%3]) >= tol:
                            isonpath = False
            
                # Includes the last point on the path
                if i == np.shape(path)[1] - 2:
                    mul = 1
                else:
                    mul = -1
            
                # Makes sure it is between the 2 points
                if any(np.array(t) > 1 + mul*tol) or any(np.array(t) < 0 - tol):
                    isonpath = False
            
                # Adds the points that passed all the tests on the path
                if isonpath:
                    onpath.append(list(p))

    newpath = unique(onpath.T) 
         
    return np.concatenate((newpath, np.ones((1, np.shape(newpath)[1]))), axis=0) #Adding ones for QE
         

def on_path_plot(path, rsc, irsc):
    epath = explicit_path(path)*np.pi*2
    epathrc = irsc.dot(epath[0:3,:])
    
    npts = path[3,:]
    path = path[:3,:]
    
    closest_points = rsc.dot(unique(np.round(epathrc))) / (2*np.pi)
    onpath = []
    dist_on_path = 0
    syms = [0] # list of positions on path of the high symmerty points
    pos_on_path = []
    for i,v in enumerate(path[:,:-1].T):
        if npts[i] > 1:
            for p in closest_points.T:
               isonpath = True
               t = []
               for j in range(3):
                   if abs((path[i+1][j]-v[j])) >= tol: # make sure there is no division by zero
                       t.append((p[j]-v[j])/(path[i+1][j]-v[j]))
                   else:
                       if abs((p[j]-v[j])) >= tol: # if so, check if nominator is near 0
                           isonpath = False
            
               # Makes sure the multiplier t is the same for each component
               if len(t) == 2:
                   if abs(t[1] - t[0]) >= tol:
                       isonpath = False
               elif len(t) == 3:
                   for j in range(3):
                       if abs(t[j]-t[(j+1)%3]) >= tol:
                           isonpath = False
            
               # Includes the last point on the path
               if i == np.shape(path)[1] - 2:
                   mul = 1
               else:
                   mul = -1
            
               # Makes sure it is between the 2 points
               if any(np.array(t) > 1 + mul*tol) or any(np.array(t) < 0 - tol):
                   isonpath = False
            
               # Adds the points that passed all the tests on the path
               if isonpath:
                   onpath.append(list(p))
                   newdist = dist_on_path + np.linalg.norm(p-v)
                   if pos_on_path == []:
                       idx = 0
                   else:
                       idx = len(pos_on_path)-1 
                       # Goes through the list to make sure it is ordered (closest point is not necessarily ordered)
                       while newdist < pos_on_path[idx]:
                           idx += -1
                   pos_on_path.insert(idx+1, dist_on_path + np.linalg.norm(p-v))

            dist_on_path += np.linalg.norm(path[i+1]-v)
            syms.append(dist_on_path)
                 
    newpath = onpath.T 
         
    return (np.concatenate((newpath, np.ones((1,np.shape(newpath)[1]))), axis=0),
            pos_on_path, syms) #Adding ones for QE

def all_points(rpc, irpc, rsc, irsc):
    # big square
    # Finds the furthest corner
    corner = np.array([[0,0,0]]).T
    for i in range(2):
        for j in range(2):
            for k in range (2):
                corner = np.reshape(np.max(np.concatenate([abs(irsc.dot(rpc).dot(np.array([[i,j,k]]).T)), corner], axis=1), axis=1),(3,1))
    
    corner = abs(np.ceil(corner)).astype(int)
    
    list_in = []
    for i in range(-corner[0,0], corner[0,0]):
        for j in range(-corner[1,0], corner[1,0]):
            for k in range (-corner[2,0], corner[2,0]):
                p = rsc.dot(np.array([i,j,k]).T)
                if (irpc.dot(p) >= 0. - tol).all() and (irpc.dot(p) < 1. - tol).all():
                    list_in.append(p.tolist())
    
    list_in = np.array(list_in).T / (2*np.pi)

    return np.concatenate((list_in, np.ones((1, np.shape(list_in)[1]))), axis=0).T #Adding ones for QE

def derivative_points(path, rsc):
    q = explicit_path(path)[:3,:] * 2 * np.pi
    
    coef = np.array([-2, -1, 1, 2])
    
    for iq, qp in enumerate(q.T):
        # For each dimention... 
        for i in range(3):
            q_disp = np.zeros((3,4))
            q_bool = [False]*4
            for k, mul in enumerate(coef):
                # ...finds the 2 nearest neighbors in each direction
                q_disp[:,k] = qp + mul*rsc[:,i]
    
            # Check if the neighbors are already in the list
            for j in range(len(q)):
                for k, mul in enumerate(coef):
                    if np.linalg.norm(q_disp[:,k] - q[:,j]) <= tol:
                        q_bool[k] = True
    
            # If none of the nieghbors is present...
            if not any(q_bool):
                # ..adds the symmetric ones
                q = np.concatenate((q.T,q_disp[:,1:3]), axis = 1)

            # If one or more neighbor is present but not next to eachother
            elif not any([a and b for a, b in zip(q_bool[:3], q_bool[1:])]):
                # finds the most centered neighbor
                q_ind = np.argmin(abs(coef[q_bool]))
                q_ind = np.where(q_bool)[0][q_ind]
                # Adds the symmetric neighbor in priority
                if q_ind < 2:
                    q = np.concatenate((q.T,[q_disp[:, q_ind+1]]), axis = 1)
                else:
                    q = np.concatenate((q.T,[q_disp[:, q_ind-1]]), axis = 1)
    
    return np.concatenate((q.T/(2*np.pi), np.ones((1, np.shape(q)[1]))), axis=0)

def to_relaxed_coord(path, irpc_perfect, rpc):
    weights = np.array(path)[3:4,:]
    path = rpc.dot(irpc_perfect).dot(path[:3,:])
    return np.concatenate((path, weights), axis=0)
