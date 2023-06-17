p5.disableFriendlyErrors = true; // disables FES

var time, dt
var nx, ny
var dx, dy

var origin

var c = 10000
var EM

function gradE(V, Avel) {
    var F = []

    for (var y = 0; y < ny; y++) {
        F[y] = []

        for (var x = 0; x < nx; x++) {

            var V00 = V[y][x];

            var V10 = (x + 1 < nx) ? V[y][x + 1] : V00;
            var V01 = (y + 1 < ny) ? V[y + 1][x] : V00;

            F[y][x] = createVector(
                - (V10 - V00) / dx - Avel[y][x].x,
                - (V01 - V00) / dy - Avel[y][x].y, 
                0
            );
        }
    }
    return F
}

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

function curl(A) {

    var F = []
    for (var y = 0; y < ny; y++) {
        F[y] = []

        for (var x = 0; x < nx; x++) {

            var A00 = A[y][x];

            var A10 = (x + 1 < nx) ? A[y][x + 1] : A00;
            var A01 = (y + 1 < ny) ? A[y + 1][x] : A00;

            F[y][x] = createVector(
                  (A01.z - A00.z) / dy,
                - (A10.z - A00.z) / dx, 
                  (A10.y - A00.y) / dx - (A01.x - A00.x) / dy
            );
        }
    }
    return F
}

class EMfield {

    constructor() {
        /* Electric field */ 
        this.V = []
        this.Vvel = []
        this.rho = []
        this.E = []
        
        /* magnetic field */
        this.A = []
        this.Avel = []
        this.J = []

        for (var i = 0; i < ny; i++) {
            /* electric field */
            this.V[i] = []
            this.Vvel[i] = []
            this.E[i] = []
            this.rho[i] = []

            /* magnetic field */
            this.A[i] = []
            this.Avel[i] = []
            this.J[i] = []

            for (var j = 0; j < nx; j++) {
                
                var x = j - gridOrigin.x 
                var y = i - gridOrigin.y

                /* e field */ 
                this.V[i][j] = 0
                this.Vvel[i][j] = 0
                this.rho[i][j] = 0
                this.E[i][j] = createVector(0, 0, 0)

                /* b field */
                this.A[i][j] = createVector(0, 0, 0)
                this.Avel[i][j] = createVector(0, 0, 0)
                this.J[i][j] = createVector(0, 0, 0)

               
            }
        }

        this.E = gradE(this.V, this.Avel)
        this.B = curl(this.A)
    }

    update(dt) {

        var acc = 0

        for (var i = 0; i < ny; i++) {
            for (var j = 0; j < nx; j++) {

                var x = j - gridOrigin.x 
                var y = i - gridOrigin.y - 30*sin(time * 1000)

                this.rho[i][j] = 10*exp(-(x*x + y*y)/0.1) 

                if (int(x) == 0 && int(y) == 0) {
                    //this.J[i][j].z = sin(time * 100)
                }
            }
        }

        for (var y = 0; y < ny; y++) {
            for (var x = 0; x < nx; x++) {

                /* electric field */ 
                var V00 = this.V[y][x];

                var V10 = (x + 1 < nx) ? this.V[y][x + 1] : V00;
                var V01 = (y + 1 < ny) ? this.V[y + 1][x] : V00;

                var V20 = (0 < x - 1) ? this.V[y][x - 1] : V00;
                var V02 = (0 < y - 1) ? this.V[y - 1][x] : V00; 
                
                acc = c * c *((V20 + V10 + V02 + V01 -4*V00) / (dx * dy) + this.rho[y][x])

                this.Vvel[y][x] += acc * dt
                this.V[y][x] += this.Vvel[y][x] * dt

                /* b field */

                var A00 = this.A[y][x];

                var A10 = (x + 1 < nx) ? this.A[y][x + 1] : A00;
                var A01 = (y + 1 < ny) ? this.A[y + 1][x] : A00;

                var A20 = (0 < x - 1) ? this.A[y][x - 1] : A00;
                var A02 = (0 < y - 1) ? this.A[y - 1][x] : A00; 
                
                /* Ax */
                acc = c * c *((A20.x + A10.x + A02.x + A01.x -4*A00.x) / (dx * dy) + this.J[y][x].x)
                this.Avel[y][x].x += acc * dt
                this.A[y][x].x += this.Avel[y][x].x * dt

                /* Ay */
                acc = c * c *((A20.y + A10.y + A02.y + A01.y -4*A00.y) / (dx * dy) + this.J[y][x].y)
                this.Avel[y][x].y += acc * dt
                this.A[y][x].y += this.Avel[y][x].y * dt

                /* Az */
                acc = c * c *((A20.z + A10.z + A02.z + A01.z -4*A00.z) / (dx * dy) + this.J[y][x].z)
                this.Avel[y][x].z += acc * dt
                this.A[y][x].z += this.Avel[y][x].z * dt
            }
        }

        this.E = gradE(this.V, this.Avel)
        this.B = curl(this.A)
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

    getB(x, y) {
        if ((x < 0) || (x > nx) || (y < 0) || (y > ny)) {
            return 0
        }
        return this.B[y][x];
    }
}

function setup() {
    createCanvas(500, 500, WEBGL)

    /* simulation utils */
    time = 0.0
    dt = 0.0

    /* sizes */
    nx = 300
    ny = floor(nx / (width/height))
    dx = width / nx
    dy = width / nx


    /* drawing utils */
    screenOrigin = createVector(width / 2, height / 2)
    gridOrigin = createVector(int(nx / 2), int(ny / 2));
    EM = new EMfield()
}

function draw() {
    dt = 0.00001
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
            let B = EM.getB(i, j)
            let f = E.mag()

            if (f > 0.1) {
                let c = color(f * 255);

                fill(c)
                
                rect(i*dx - screenOrigin.x, j*dy - screenOrigin.y, dx, dy)
            }
            
            
        }
    }
}