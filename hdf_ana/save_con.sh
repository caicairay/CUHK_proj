DIRECTORY=../figs/E_avr
SWITCHER=3 #1=[den & b] ; 2=[Ek & Eb, contourf]; 3=[Ek & Eb, avr]; 4=[vr & vt, avr] 5=[E-Ek-Eb]
GIFNAME=E_avr
up_lim=300
low_lim=0
FILE=scp_hdf.py
LOGFILE="$FILE".log
if [ -d "$DIRECTORY" ]&&[ "$low_lim" -eq 0 ]; then
	if test "$(ls -A "$DIRECTORY")"; then
		echo directory not empty
	else
		for ((c=low_lim;c<=up_lim;c++))
		do
			python -W ignore "$FILE" $c $SWITCHER "$DIRECTORY" >> "$LOGFILE" 
		done
                convert -resize 40% -delay 10 -loop 0 "$DIRECTORY"/*.png "$DIRECTORY"/"$GIFNAME".gif
                scp "$DIRECTORY"/"$GIFNAME".gif zcao@zcaooffice.hopto.org:~/Desktop
	fi
else
	mkdir -p "$DIRECTORY"
	for ((c=low_lim;c<=up_lim;c++))
	do
		python -W ignore "$FILE" $c $SWITCHER "$DIRECTORY" >> "$LOGFILE"
	done
        convert -resize 40% -delay 10 -loop 0 "$DIRECTORY"/*.png "$DIRECTORY"/"$GIFNAME".gif
        scp "$DIRECTORY"/"$GIFNAME".gif zcao@zcaooffice.hopto.org:~/Desktop
fi

