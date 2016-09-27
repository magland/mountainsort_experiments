// Remove all rejected tags
for (var i in clusters) {
	clusters[i].removeTag('rejected');
}

console.log('peak_amp');
for (var i in clusters) {	
	if (clusters[i].metric('peak_amp')<5)
		clusters[i].addTag('rejected');
}

console.log('peak_noise');
for (var i in clusters) {	
	if (clusters[i].metric('peak_noise')>30)
		clusters[i].addTag('rejected');
}

console.log('firing_rate');
for (var i in clusters) {	
	if (clusters[i].metric('firing_rate')<0.1)
		clusters[i].addTag('rejected');
}
