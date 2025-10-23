.PHONY: matrix m2s1_tx
	

matrix:
	mkdir -p ~/mnt/gmt-full-field-display/matrix/
	for id in {1..7}; do cargo r -r -p gmt-full-field-display -- --n-thread 25 segment --id $$id --tx=100; done
	mv gmt-full-field_M2S* ~/mnt/gmt-full-field-display/matrix/
	for id in {1..7}; do cargo r -r -p gmt-full-field-display -- --n-thread 25 segment --id $$id --ty=100; done
	mv gmt-full-field_M2S* ~/mnt/gmt-full-field-display/matrix/
	for id in {1..7}; do cargo r -r -p gmt-full-field-display -- --n-thread 25 segment --id $$id --rx=10; done
	mv gmt-full-field_M2S* ~/mnt/gmt-full-field-display/matrix/
	for id in {1..7}; do cargo r -r -p gmt-full-field-display -- --n-thread 25 segment --id $$id --ry=-10; done
	mv gmt-full-field_M2S* ~/mnt/gmt-full-field-display/matrix/

m2s1_tx:
	mkdir -p ~/mnt/gmt-full-field-display/M2S1_tx
	for tx in {-100..100..20}; do cargo r -r -p gmt-full-field-display -- --n-thread 25 -f M2S1_tx$$tx.pkl segment --id 1 --tx=$$tx; done
	mv M2S1_tx* ~/mnt/gmt-full-field-display/M2S1_tx/
