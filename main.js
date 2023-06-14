
var time, dt
var nx, ny
var dx, dy

var origin

var EM

function gradient(V) {

    var F = []

    for (var y = 0; y < ny; y++) {
        F[y] = []

        for (var x = 0; x < nx; x++) {

            var V00 = V[y][x];

            var V10 = (x + 1 < nx) ? V[y][x + 1] : V00;
            var V01 = (y + 1 < ny) ? V[y + 1][x] : V00;

            F[y][x] = createVector(
                - (V10 - V00) / dx,
                - (V01 - V00) / dy, 
                0
            );
        }
    }
    return F
}

class EMfield {

    constructor() {
        this.V = []
        this.Vvel = []

        this.rho = []

        this.E = []

        for (var i = 0; i < ny; i++) {
            this.V[i] = []
            this.Vvel[i] = []

            this.rho[i] = []
            
            this.E[i] = []
            for (var j = 0; j < nx; j++) {
                
                var x = j - gridOrigin.x 
                var y = i - gridOrigin.y

                this.V[i][j] = 0
                this.Vvel[i][j] = 0

                this.rho[i][j] = 0

                this.E[i][j] = createVector(0, 0, 0)
            }
        }

        this.E = gradient(this.V)

    }

    update(dt) {

        var acc = 0, vel = 0

        for (var i = 0; i < ny; i++) {
            for (var j = 0; j < nx; j++) {

                var x = j - gridOrigin.x - sin(time)
                var y = i - gridOrigin.y

                this.rho[i][j] = exp(-(x*x + y*y)/0.1)
            }
        }

        for (var y = 0; y < ny; y++) {
            vel = 0
            for (var x = 0; x < nx; x++) {

                var V00 = this.V[y][x];

                var V10 = (x + 1 < nx) ? this.V[y][x + 1] : V00;
                var V01 = (y + 1 < ny) ? this.V[y + 1][x] : V00;

                var V20 = (0 < x - 1) ? this.V[y][x - 1] : V00;
                var V02 = (0 < y - 1) ? this.V[y - 1][x] : V00; 
                
                acc = (V20 + V10 + V02 + V01 -4*V00) / (dx * dy) + this.rho[y][x]

                this.Vvel[y][x] += acc * dt
                this.V[y][x] += this.Vvel[y][x] * dt
            }
        }

        this.E = gradient(this.V)
    }

    getE(x, y) {
        if ((x < 0) || (x > nx) || (y < 0) || (y > ny)) {
            return createVector(0, 0, 0)
        }
        return this.E[y][x];
    }

    getV(x, y) {
        if ((x < 0) || (x > nx) || (y < 0) || (y > ny)) {
            return 0
        }
        return this.V[y][x];
    }

    setRho(x, y, value) {
        this.rho[y][x] = value
    }
}

function setup() {
    createCanvas(windowWidth, windowHeight, WEBGL)

    /* simulation utils */
    time = 0.0
    dt = 0.0

    /* sizes */
    nx = 70
    ny = floor(nx / (width/height))
    dx = width / nx
    dy = width / nx


    /* drawing utils */
    screenOrigin = createVector(width / 2, height / 2)
    gridOrigin = createVector(int(nx / 2), int(ny / 2));
    EM = new EMfield()
}

function draw() {
    dt = 0.1
    time += dt
    clear()
    background(0)
    
    EM.update(dt)

    /* draw */
    for (var i = 0; i < nx; i++) {
        for (var j = 0; j < ny; j++) {
            var x = i - gridOrigin.x
            var y = j - gridOrigin.y
            

            let E = EM.getE(i, j)
            let f = E.mag()
            let c = color(f * 255);

            fill(c)

            rect(i*dx - screenOrigin.x, j*dy - screenOrigin.y, dx, dy)
            
        }
    }
}