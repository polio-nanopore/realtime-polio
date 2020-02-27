
## Running RAMPART

With this set up, to run RAMPART:

```
rampart --protocol path/to/realtime-polio/rampart 
```

where path/to/realtime-polio/rampart will change depending on each machine and where you have saved the realtime-polio directory.

If you navigate to realtime-polio/rampart in your computer, then type 

``pwd``

this will tell you the path to the directory. Copy and paste this into your terminal when you're running RAMPART. 

> *Note*: This command will not change in between runs. Once you know what your path is, you can use the same command every time. For instance, on my machine, I type: ``rampart --protocol /Users/aine/repositories/realtime-polio/rampart``. Yours will certainly be something similar. The difference in between runs is *where* you are when you type this command. You must be in the directory you have set up for the current MinION run. 

Open a web browser to view [http://localhost:3000](http://localhost:3000)

More information about RAMPART can be found [here](https://github.com/artic-network/rampart).

## RAMPART command line options

```
usage: rampart [-h] [-v] [--verbose] [--ports PORTS PORTS]
               [--protocol PROTOCOL] [--title TITLE]
               [--basecalledPath BASECALLEDPATH]
               [--annotatedPath ANNOTATEDPATH]
               [--referencesPath REFERENCESPATH]
               [--referencesLabel REFERENCESLABEL]
               [--barcodeNames BARCODENAMES [BARCODENAMES ...]]
               [--annotationOptions ANNOTATIONOPTIONS [ANNOTATIONOPTIONS ...]]
               [--clearAnnotated] [--simulateRealTime SIMULATEREALTIME]
               [--devClient] [--mockFailures]
```

## Stopping RAMPART

If you need to stop RAMPART, you can type ``Ctrl+c`` in the terminal window running RAMPART.
