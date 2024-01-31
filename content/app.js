// const server = "http://tebas.lpi.tel.uva.es:8080"
const server = "http://127.0.0.1:8080"

function callKomaMRI(){
    var phantom = document.getElementById("phantom");
    canvas.height = 0; canvas.width = 0;

    // Input parameters: mat (-> Sequence) & vec (-> Scanner)
    // This parameters will come from the sequence editor (GUI) in future versions
    var mat  = [[1,        2,      5],           // cod
                [5.87e-4,  0.01,   0],           // dur
                [0,        0,      0],           // gx
                [0,        0,      0],           // gy
                [1,        0,      0],           // gz
                [10e-6,    0,      0],           // b1x
                [0,        0,      0],           // b1y
                [0,        0,      0],           // Δf
                [0,        0,      0.4],         // fov
                [0,        0,      201]];        // n

    var vec  =  [1.5,          //B0
                 10e-6,        //B1
                 2e-6,         //Delta_t
                 60e-3,        //Gmax
                 500];         // Smax

    var params = {
        phantom: phantom.value,
        mat: mat,
        vec: vec
    }

    // HTTP Status Codes:
    // 200: OK
    // 202: Accepted
    // 303: See other

    document.getElementById("simButton").disabled = true;

    fetch(server + "/simulate",{
        method: "POST",
        headers:{
            "Content-type": "application/json",
        },
        body: JSON.stringify(params)})
    .then(res => {
            if ((res.status == 202) && (loc = res.headers.get('location'))){
                requestResult(loc)
            }else{
                // Error
            }
        }
    )
}


function requestResult(loc){
    fetch(server + loc)
        .then(res => {
            if(res.redirected){
            // Redirected indica que la respuesta es fruto de una redirección.
            // Esto quiere decir que el cliente ha recibido un 303 y ha hecho
            // automáticamente la petición a simulate/{simulationId}/status
                document.getElementById('response').style.visibility = "visible";
                res.json()
                .then(json => {
                    if(json == -1){
                        document.getElementById('response').innerHTML = "Starting simulation...";
                    } else {
                        // Status Bar
                        document.getElementById('response').innerHTML = json + "%";

                        document.getElementById('myProgress').style.visibility = "visible";
                        var elem = document.getElementById("myBar");
                        elem.style.width = json + "%";
                    }
                    if(json<100){
                        setTimeout(requestResult(loc),500);
                    } else if (json == 100) {
                        // Status Bar
                        document.getElementById('response').innerHTML = "Reconstructing...";
                        document.getElementById('myProgress').style.visibility = "collapse";
                        setTimeout(requestResult(loc),500);
                    }
                })
            }else{
            // Aquí estaríamos obteniendo el resultado de la simulación
            // El cliente ha recibido un 200 como respuesta al GET simulate/{simulationId}
                let image = []
                res.json()
                .then(json => json.data)
                .then(function(array){
                    N = Math.sqrt(array.length); // Suponemos imágenes cuadradas

                    var canvas = document.getElementById( 'canvas' );
                    canvas.height = N; canvas.width = N;
                    var context = canvas.getContext( '2d' );

                    var imgArray = []
                    for (var i = 0; i < array.length; i++){
                        imgArray[4*i] = array[i];
                        imgArray[4*i+1] = array[i];
                        imgArray[4*i+2] = array[i];
                        imgArray[4*i+3] = 255;
                    }

                    document.getElementById('response').style.visibility = "collapse";

                    document.getElementById('myProgress').style.visibility = "collapse";
                    document.getElementById("myBar").style.width = 0 + "%";

                    const imgData = new ImageData(Uint8ClampedArray.from(imgArray), N, N);
                    context.putImageData(imgData, 0, 0);

                    document.getElementById("simButton").disabled = false;
                })
            }
        })
}
