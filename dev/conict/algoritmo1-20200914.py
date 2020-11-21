        t=t.tolist()
        v=v.tolist()
        s=s.tolist()
        f=f.tolist()
        a=a.tolist()
        Em=Em.tolist()
        
        t_p=t_p.tolist()
        v_p=v_p.tolist()
        s_p=s_p.tolist()
        f_p=f_p.tolist()
        a_p=a_p.tolist()
        Em_p=Em_p.tolist()
        
        if len(t)-1 != len(t_p)-1:
            n=max((len(t)-1), len(t_p)-1)
            if len(t)-1 > len(t_p)-1:
                for i in range(n-(len(t_p)-1)):
                    
                    t_p.append(nan)
                    
                    v_p.append(nan)
                    
                    s_p.append(nan)
                    
                    f_p.append(nan)
                    
                    a_p.append(nan)
                    
                    Em_p.append([nan, nan, nan, nan, nan])
                    
            if len(t)-1 < len(t_p)-1:
                for j in range(n-(len(t)-1)):
                    
                    t.append(nan)
                    
                    v.append(nan)
                    
                    s.append(nan)
                    
                    f.append(nan)
                    
                    a.append(nan)
                    
                    Em.append([nan, nan, nan, nan, nan])
                    
        t=array(t)
        v=array(v)
        s=array(s)
        f=array(f)
        a=array(a)
        Em=array(Em)
        
        t_p=array(t_p)
        v_p=array(v_p)
        s_p=array(s_p)
        f_p=array(f_p)
        a_p=array(a_p)
        Em_p=array(Em_p)