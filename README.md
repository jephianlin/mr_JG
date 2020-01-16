**This code is a copy of  Minimum Rank Library from Jason Grout's GitHub. 

Minimum Rank Sage Library
=========================

One goal of this library is to be able to be used both by loading the necessary files in the Sage notebook via load statements, as well as installation as a python module in Sage.  Thus, it is important that the names in the separate submodules do not trample on each other (since the load command imports everything at the top level)

See http://sage.cs.drake.edu/home/pub/69/ for an example of how to load and use this library in Sage.  In particular, the following code in a Sage notebook cell will load this library:

```python
URL='https://raw.githubusercontent.com/jephianlin/mr_JG/master/'
### *.pyx files failed to compile in newer Sage
### Will take some time to fix it
# files=['Zq_c.pyx','Zq.py','zero_forcing_64.pyx','zero_forcing_wavefront.pyx','minrank.py', 'inertia.py']
files=['Zq.py','minrank.py', 'inertia.py']
for f in files:
    print("Loading %s..."%f);
    load(URL+f)
```
  
CoCalc user
-----------

For CoCalc (previously SageMathCloud) user, if you are using a free acount, then the code above does not work for you.  It is because the connection from CoCalc to outside is forbidden.  You may do the following:

1. Click the green button "Clone or download" and choose "Download Zip".
2. Go to your project in CoCalc and upload this zip file.
3. Create a terminal and open it.
4. Type in the terminal::
```bash
unzip mr_JG-master.zip
```
5. Open your Sage worksheet and execute::
```python
URL='mr_JG-master/'
### *.pyx files failed to compile in newer Sage
### Will take some time to fix it
# files=['Zq_c.pyx','Zq.py','zero_forcing_64.pyx','zero_forcing_wavefront.pyx','minrank.py', 'inertia.py']
files=['Zq.py','minrank.py', 'inertia.py']
for f in files:
    print("Loading %s..."%f);
    load(URL+f)
```
If you are tired of downloading and uploading all the time, you may consider to subscribe CoCalc, which give you internet access from CoCalc to outside (in particular, GitHub).
