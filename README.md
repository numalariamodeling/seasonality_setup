# seasonality_setup
combining Ben's bf-seasonal fitting code with tobias' project-template 

**Workflow:**  
1. Run param_select.R with rnd=0  
2. Run fitting.py  woth rnd=0  
3. Run param_select.R with rnd=1  
4. Run fitting.py  woth rnd=1  
5. Keep running with rnd=rnd+1 until satisfied.
