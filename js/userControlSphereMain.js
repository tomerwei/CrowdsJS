import * as THREE from "three";
import * as PHY from "simplePhysics";
import {OrbitControls} from "three/addons/controls/OrbitControls.js";

import Stats from "three/addons/libs/stats.module.js";

let renderer, scene, camera;
let world = {
    x: 80,
    z: 80,
};
let agentData = [];
let pickableObjects = [];
let selected = [];
let mouse = new THREE.Vector2();
const raycaster = new THREE.Raycaster();
let grid, ring;

let spotLights = {};
let topTextures = {};
let topTexture;
const RADIUS = 1;
const blueAgentMaterial = new THREE.MeshLambertMaterial({
    color: 0x0000ff,
});
const redAgentMaterial = new THREE.MeshLambertMaterial({
    color: 0xff0000,
});
const greenAgentMaterial = new THREE.MeshLambertMaterial({
    color: 0x00ff00,
});
const yellowAgentMaterial = new THREE.MeshLambertMaterial({
    color: 0xffff00,
});

const stats = new Stats();
document.body.appendChild(stats.dom);

let arrows = [];
let g_arrows = [];
let xarrows = [];
let zarrows = [];

let parameters = {
    best:[],
    wallData: [],
}

const WORLDUNIT = 1
const tile = {
    w:WORLDUNIT * 2,
    h:WORLDUNIT * 2
}
let tiles = [];
let wallData = [];
class Tile {

    constructor( r, c, x, y, z, cost = 1, weight = 0) {

        this.x = x;
        this.y = y;
        this.z = z;

        this.r = r;
        this.c = c;

        this.cost = cost;

        this.g = cost;
        this.h = 0;
        this.f = 0;

        this.parent = null;

        this.weight = weight
        this.FFNeighbors = []
        this.FFMoreNeighbors = []
        this.vec = {x:0, z:0}

        this.density = 0;
    }
}

function cut(){
    let world_width = world.x;
    let world_height = world.z;

    let tile_width = tile.w;
    let tile_height = tile.h;

    let Rs =   Math.floor(world_width / tile_width);
    let Cs = Math.floor(world_height / tile_height);

    if (world_width % tile_height !== 0){
        Rs += 1;
    }

    if (world_height % tile_height !== 0){
        Cs += 1;
    }

    return [Rs, Cs];
}

function gridization(){

    [rows, columns] = cut();


    const start_point = {
        x: 0 - world.x / 2,
        y: 0,
        z: 0 - world.z / 2,
    };



    for (let i = 0; i < rows; i++) {
        tiles[i] = [];



        for (let j = 0; j < columns; j++) {

            const object_position = {
                x: start_point.x + WORLDUNIT + i * tile.w,
                y: 1,
                z: start_point.z + WORLDUNIT + j * tile.h,
            };

            let cost;
            if (checkContainsTuple(obstacles, [i, j])){
                cost = obstacleCost;
            }else {
                cost = normalCost;

            }


            tiles[i][j] = new Tile(i, j, object_position.x, object_position.y, object_position.z, cost);


        }
    }

    // // set flow field
    // let columnWidth = tiles[0].length;
    // let vectors = makeFlowField();
    // // console.log(vectors);
    // for (const row of vectors) {
    //     for (const column of row){
    //         let cell = column;
    //         let x = cell.threeVec[0];
    //         let z = cell.threeVec[1];
    //         tiles[cell.r][cell.c].vec = {x:x, z:z};
    //         tiles[cell.r][cell.c].digit = cell.digit;
    //     }
    // }

    // console.log(tiles);

    // for (let i = 0; i< rows;i++){
    //     for (let j = 0; j<columns;j++){
    //         // Create a box geometry and a mesh material
    //
    //
    //         let t = tiles[i][j];
    //
    //         const geometry = new THREE.BoxGeometry(tile.w, WORLDUNIT * 2, tile.h);
    //
    //         let material;
    //         if (t.cost >= obstacleCost){
    //             material = new THREE.MeshStandardMaterial({
    //                 transparent: true,
    //                 opacity: 1.0,
    //                 color: 0x333333 // set a color to disable the transparent effect
    //             });
    //
    //
    //         }else {
    //             material = new THREE.MeshStandardMaterial({
    //                 transparent: true,
    //                 opacity: 0.0,
    //                 color: 0x00ff00 // set a color to see the transparent effect
    //             });
    //         }
    //
    //
    //         // Create a mesh by combining the geometry and the material
    //         const cube = new THREE.Mesh(geometry, material);
    //
    //         // Set the mesh's name and userData properties
    //         cube.name = "MyCube_" + i + "_" + j;
    //         cube.userData = {
    //             type: "box",
    //             x: t.x,
    //             y: t.y,
    //             z: t.z,
    //             r: t.r,
    //             c: t.c,
    //         };
    //         cube.position.set(t.x, t.y, t.z);
    //
    //
    //         if(t.cost >= obstacleCost){
    //             pickableWall.push(cube);
    //             wallData.push(cube.userData);
    //
    //         }else {
    //             pickableWalkingTiles.push(cube);
    //         }
    //
    //
    //         pickableTiles.push(cube);
    //
    //         // Add the mesh to the scene
    //         scene.add(cube);
    //
    //
    //         // break;
    //     }
    //     // break;
    // }



    // add exit grids
    // opens.forEach(function (e, index){
    //
    //     exits[index].x = tiles[e[0]][e[1]].x
    //     exits[index].z = tiles[e[0]][e[1]].z
    //
    // });
    //
}


function init() {
    // renderer
    renderer = new THREE.WebGLRenderer();
    renderer.shadowMap.enabled = true;
    renderer.shadowMap.type = THREE.PCFSoftShadowMap; //
    renderer.setSize(window.innerWidth, window.innerHeight);
    document.body.appendChild(renderer.domElement);

    // scene
    scene = new THREE.Scene();
    // camera
    camera = new THREE.PerspectiveCamera(
        45,
        window.innerWidth / window.innerHeight,
        1,
        1000
    );

    camera.position.set(-67.26, 54.07, -3.77);
    camera.rotation.order = "YXZ";
    camera.rotation.y = -1.6267;
    camera.rotation.x = -0.46;

    // controls
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.addEventListener("change", render);
    controls.enableZoom = false;
    controls.enablePan = false;
    controls.maxPolarAngle = Math.PI / 2;

    // light
    const light = new THREE.PointLight(0xffffff, 0.9, 0, 100000);
    light.position.set(0, 50, 120);
    light.castShadow = true;
    light.shadow.mapSize.width = 512; // default
    light.shadow.mapSize.height = 512; // default
    light.shadow.camera.near = 0.5; // default
    light.shadow.camera.far = 5000; // default

    const directionalLight = new THREE.DirectionalLight(0xffffff, 0.5);
    directionalLight.castShadow = true;
    directionalLight.position.set(-5, 20, 4);
    directionalLight.target.position.set(9, 0, -9);
    directionalLight.shadow.camera.left *= 9;
    directionalLight.shadow.camera.right *= 9;
    directionalLight.shadow.camera.top *= 9;
    directionalLight.shadow.camera.bottom *= 9;

    scene.add(directionalLight);

    // axes
    scene.add(new THREE.AxesHelper(40));
    const loader = new THREE.TextureLoader();
    const texture = loader.load("resources/OIP.jpg");
    texture.wrapS = THREE.RepeatWrapping;
    texture.wrapT = THREE.RepeatWrapping;
    texture.magFilter = THREE.NearestFilter;
    const repeats = 40 / 32;
    texture.repeat.set(repeats, repeats);

    topTexture = loader.load("resources/triangle2.png");
    //topTexture.wrapS = THREE.RepeatWrapping;
    //topTexture.wrapT = THREE.RepeatWrapping;
    topTexture.magFilter = THREE.NearestFilter;
    topTexture.repeat.set(3, 3);
    //topTexture.rotation = -Math.PI / 2;
    // grid
    const geometry = new THREE.PlaneGeometry(world.x, world.z, 10, 10);
    const material = new THREE.MeshPhongMaterial({
        map: texture,
        //side: THREE.DoubleSide,
    });
    grid = new THREE.Mesh(geometry, material);
    grid.castShadow = true; //default is false
    grid.receiveShadow = true; //default
    grid.rotation.order = "YXZ";
    grid.rotation.y = -Math.PI / 2;
    grid.rotation.x = -Math.PI / 2;
    scene.add(grid);

    // experiment border
    // const boxGeometry1 = new THREE.BoxGeometry(50, 5, 1);
    // const boxMaterial1 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
    // const left = new THREE.Mesh(boxGeometry1, boxMaterial1);
    // left.position.set(0, 2.5, -25);
    //
    // const boxGeometry2 = new THREE.BoxGeometry(1, 5, 50);
    // const boxMaterial2 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
    // const bottom = new THREE.Mesh(boxGeometry2, boxMaterial2);
    // bottom.position.set(-25, 2.5, 0);
    //
    // const boxGeometry3 = new THREE.BoxGeometry(1, 5, 50);
    // const boxMaterial3 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
    // const top = new THREE.Mesh(boxGeometry3, boxMaterial3);
    // top.position.set(25, 2.5, 0);
    //
    // const boxGeometry4 = new THREE.BoxGeometry(50, 5, 1);
    // const boxMaterial4 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
    // const right = new THREE.Mesh(boxGeometry4, boxMaterial4);
    // right.position.set(0, 2.5, 25);
    //
    // scene.add(left);
    // scene.add(bottom);
    // scene.add(top);
    // scene.add(right);



    const ringGeometry = new THREE.RingGeometry(1, 3, 12);
    const ringMaterial = new THREE.MeshBasicMaterial({
        color: 0xffff00,
        side: THREE.DoubleSide,
    });
    ring = new THREE.Mesh(ringGeometry, ringMaterial);
    scene.add(ring);
    ring.rotation.x = -Math.PI / 2;
    ring.position.y += 0.01;


    function sampleCirclePoints(radius, sampleCount, centerX, centerY) {
        let points = [];
        for (let i = 0; i < sampleCount; i++) {
            // Angle in radians
            let angle = 2 * Math.PI * i / sampleCount;

            // Calculating x and y coordinates
            let x = centerX + radius * Math.cos(angle);
            let y = centerY + radius * Math.sin(angle);

            points.push({ x: x, y: y });
        }
        return points;
    }


    function sampleCirclePointsWithDistance(sampleCount, distanceBetweenPoints, centerX, centerY) {
        // Calculating the total circumference that would fit the points with the given distance
        let totalCircumference = distanceBetweenPoints * sampleCount;

        // Updating the radius based on the new circumference
        let radius = totalCircumference / (2 * Math.PI);

        let points = [];
        for (let i = 0; i < sampleCount; i++) {
            // Angle in radians
            let angle = 2 * Math.PI * i / sampleCount;

            // Calculating x and y coordinates
            let x = centerX + radius * Math.cos(angle);
            let y = centerY + radius * Math.sin(angle);

            points.push({ x: x, y: y });
        }
        return points;
    }



    function testScenario(){
        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 4,
            {
                x: 20,
                z: 0,
            },
            {
                x: -20,
                z: -5,
            },
            0.8,
            "X"
        );

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 4,
            {
                x: -20,
                z: 0,
            },
            {
                x: 20,
                z: 0,
            },
            0.8,
            "X"
        );
    }

    function testHallwayScenario(){

        for(let i=0;i<1;i++){
            addColumnAgentGroup(
                agentData,
                1,
                RADIUS * 1.5,
                {
                    x: 20 ,
                    z: 0 + i * 6,
                },
                {
                    x: -20,
                    z: 0 + i * 6,
                },
                0.8,
                "X"
            );

            addColumnAgentGroup(
                agentData,
                1,
                RADIUS * 4,
                {
                    x: -20,
                    z: 0 + i * 6 - 2
                },
                {
                    x: 20,
                    z: 0 + i * 6  - 2
                },
                0.8,
                "X"
            );


        }


    }

    function testCrossScenario(){
        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 4,
            {
                x: 20,
                z: 0,
            },
            {
                x: -20,
                z: -5,
            },
            0.8,
            "X"
        );

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 4,
            {
                x: 0,
                z: 20,
            },
            {
                x: 0,
                z: -20,
            },
            0.8,
            "X"
        );
    }
    function testCrossWithDiagnoScenario(){
        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 4,
            {
                x: 20,
                z: 0,
            },
            {
                x: -20,
                z: -5,
            },
            0.8,
            "X"
        );

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 4,
            {
                x: -10,
                z: 10,
            },
            {
                x: 20,
                z: -20,
            },
            0.8,
            "X"
        );
    }

    function defaultScenario(){
        addColumnAgentGroup(
            agentData,
            50,
            RADIUS * 4,
            {
                x: 0,
                z: 0,
            },
            {
                x: -10,
                z: 10,
            },
            0.8,
            "X"
        );

    }

    function circleScenario(){

        let points = sampleCirclePoints(20, 20, 0, 0);
        // let points = sampleCirclePointsWithDistance(42, 2 * 2* RADIUS + 2, 0, 0);
        console.log(points);
        points.forEach(function (point){
            addColumnAgentGroup(
                agentData,
                1,
                0,
                {
                    x: point.x,
                    z: point.y,
                },
                {
                    x: -point.x ,
                    z: -point.y,
                },
                0.8,
                "X"
            );
        });

    }

    //********************************************************* Trying scenario 1 *********************************************************
    // Time: 0.06 seconds of the video

    function tryingScenario_Bilas_1() {
        for (let i = 0; i < 5; i++) {
            for (let j = 0; j < 3; j++) {
                addColumnAgentGroup(
                    agentData,
                    1,
                    RADIUS * 1.5,
                    {
                        x: 30 - i * 6,
                        //x: 30,
                        z: -10 + j * 6,
                    },
                    {
                        x: -999,
                        //x: -20,
                        z: -10 + j * 6,
                    },
                    0.8,
                    "X"
                );

                addColumnAgentGroup(
                    agentData,
                    1,
                    RADIUS * 1.5,
                    {
                        x: -30 + i * 6,
                        //x: 30,
                        z: -6 + j * 6,
                    },
                    {
                        x: 999,
                        //x: -20,
                        z: -6 + j * 6,
                    },
                    0.8,
                    "X"
                );

            }
        }
    }

    function tryingScenario_Bilas_1_debug() {

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 1.5,
            {
                x: 30,
                z: -10,
            },
            {
                x: -999,
                z: -10,
            },
            0.8,
            "X"
        );

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 1.5,
            {
                x: -30,
                //x: 30,
                z: -6,
            },
            {
                x: 999,
                z: -6,
            },
            0.8,
            "X"
        );


        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 1.5,
            {
                x: 30,
                z: -10 + 1 * 6,
            },
            {
                x: -999,
                z: -10 + 1 * 6,
            },
            0.8,
            "X"
        );

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 1.5,
            {
                x: -30,
                //x: 30,
                z: -6 + 1 * 6,
            },
            {
                x: 999,
                //x: -20,
                z: -6 + 1 * 6,
            },
            0.8,
            "X"
        );


        // for (let i = 0; i < 1; i++) {
        //     for (let j = 0; j < 2; j++) {
        //         addColumnAgentGroup(
        //             agentData,
        //             1,
        //             RADIUS * 1.5,
        //             {
        //                 x: 30 - i * 6,
        //                 //x: 30,
        //                 z: -10 + j * 6,
        //             },
        //             {
        //                 x: -999,
        //                 //x: -20,
        //                 z: -10 + j * 6,
        //             },
        //             0.8,
        //             "X"
        //         );
        //
        //         addColumnAgentGroup(
        //             agentData,
        //             1,
        //             RADIUS * 1.5,
        //             {
        //                 x: -30 + i * 6,
        //                 //x: 30,
        //                 z: -6 + j * 6,
        //             },
        //             {
        //                 x: 999,
        //                 //x: -20,
        //                 z: -6 + j * 6,
        //             },
        //             0.8,
        //             "X"
        //         );
        //
        //     }
        // }
    }

    function tryingScenario_sphere_Bilas_1_debug() {

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 1.5,
            {
                x: 30,
                z: -11,
            },
            {
                x: -999,
                z: -11,
            },
            0.8,
            "X"
        );

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 1.5,
            {
                x: -30,
                z: -7,
            },
            {
                x: 999,
                z: -7,
            },
            0.8,
            "X"
        );


        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 1.5,
            {
                x: 30,
                z: -7,
            },
            {
                x: -999,
                z: -7,
            },
            0.8,
            "X"
        );

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 1.5,
            {
                x: -30,
                //x: 30,
                z: -3,
            },
            {
                x: 999,
                //x: -20,
                z: -3,
            },
            0.8,
            "X"
        );


        // for (let i = 0; i < 1; i++) {
        //     for (let j = 0; j < 2; j++) {
        //         addColumnAgentGroup(
        //             agentData,
        //             1,
        //             RADIUS * 1.5,
        //             {
        //                 x: 30 - i * 6,
        //                 //x: 30,
        //                 z: -10 + j * 6,
        //             },
        //             {
        //                 x: -999,
        //                 //x: -20,
        //                 z: -10 + j * 6,
        //             },
        //             0.8,
        //             "X"
        //         );
        //
        //         addColumnAgentGroup(
        //             agentData,
        //             1,
        //             RADIUS * 1.5,
        //             {
        //                 x: -30 + i * 6,
        //                 //x: 30,
        //                 z: -6 + j * 6,
        //             },
        //             {
        //                 x: 999,
        //                 //x: -20,
        //                 z: -6 + j * 6,
        //             },
        //             0.8,
        //             "X"
        //         );
        //
        //     }
        // }
    }

//********************************************************* Trying scenario 1 *********************************************************


//********************************************************* Trying scenario 2 *********************************************************
    // Time: 3.06 seconds seconds of the video
    function tryingScenario_Bilas_2() {

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 1.5,
            {
                x: 20,
                z: 0,
            },
            {
                x: -20,
                z: 0,
            },
            0.8,
            "X"
        );
        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 1.5,
            {
                x: 20,
                z: 0,
            },
            {
                x: -20,
                z: 0,
            },
            0.8,
            "X"
        );
        //-----------------------------------------------

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 4,
            {
                x: -20,
                z: 0
            },
            {
                x: 20,
                z: 0
            },
            0.8,
            "X"
        );

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 4,
            {
                x: -20,
                z: 0
            },
            {
                x: 20,
                z: 0
            },
            0.8,
            "X"
        );
    }
//********************************************************* Trying scenario 2 *********************************************************

//********************************************************* Trying scenario 3 *********************************************************
    // Time: 5:11 seconds of the video: https://youtu.be/iC8KHkoZR8k?t=314

    function tryingScenario_Bilas_3() {

        for (let i = 0; i < 5; i++) {
            for (let j = 0; j < 3; j++) {
                addColumnAgentGroup(
                    agentData,
                    1,
                    RADIUS * 1.5,
                    {
                        x: 10 - i * 6,
                        //x: 30,
                        z: 20 - j * 6,
                    },
                    {
                        x: 10 - i * 6,
                        //x: -20,
                        z: -999,
                    },
                    0.8,
                    "X"
                );

            }
        }

        // experiment border
        const boxGeometry1 = new THREE.BoxGeometry(31, 5, 1);
        const boxMaterial1 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
        const left = new THREE.Mesh(boxGeometry1, boxMaterial1);
        left.position.set(20, 2.5, 0);

        wallData.push({
            center: new THREE.Vector3(20, 2.5, 0),
            depth: 1, // along z-axis
            width: 31, // along x-axis
            base: new THREE.Vector3(20 - 31/2, 2.5, 0),
            tip: new THREE.Vector3(20 + 31/2, 2.5, 0),
        });

        const boxGeometry2 = new THREE.BoxGeometry(1, 5, 50);
        const boxMaterial2 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
        const bottom = new THREE.Mesh(boxGeometry2, boxMaterial2);
        //bottom.position.set(-15, 2.5, 0);
        bottom.position.set(-15, 2.5, -25);

        wallData.push({
            center: new THREE.Vector3(-15, 2.5, -25),
            depth: 50, // along z-axis
            width: 1, // along x-axis
            base: new THREE.Vector3(-15, 2.5, -25 - 50/2),
            tip: new THREE.Vector3(-15, 2.5, -25 + 50 / 2),
        });

        const boxGeometry3 = new THREE.BoxGeometry(1, 5, 50);
        const boxMaterial3 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
        const top = new THREE.Mesh(boxGeometry3, boxMaterial3);
        top.position.set(5, 2.5, -25);
        wallData.push({
            center: new THREE.Vector3(5, 2.5, -25),
            depth: 50, // along z-axis
            width: 1, // along x-axis
            base: new THREE.Vector3(5, 2.5, -25 - 50/2),
            tip: new THREE.Vector3(5, 2.5, -25 + 50/2),
        });

        const boxGeometry4 = new THREE.BoxGeometry(31, 5, 1);
        const boxMaterial4 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
        const right = new THREE.Mesh(boxGeometry4, boxMaterial4);
        right.position.set(-30, 2.5, 0);
        wallData.push({
            center: new THREE.Vector3(-30, 2.5, 0),
            depth: 1, // along z-axis
            width: 31, // along x-axis
            base: new THREE.Vector3(-30 - 31 / 2, 2.5, 0),
            tip: new THREE.Vector3(-30 + 31 / 2, 2.5, 0),
        });

        scene.add(left);
        scene.add(bottom);
        scene.add(top);
        scene.add(right);

        parameters.wallData = wallData;

    }

    function tryingScenario_Bilas_3_debug() {

        // addColumnAgentGroup(
        //     agentData,
        //     1,
        //     RADIUS * 1.5,
        //     {
        //         x: 9,
        //         //x: 30,
        //         z: 14,
        //     },
        //     {
        //         x: 9,
        //         //x: -20,
        //         z: -999,
        //     },
        //     0.8,
        //     "X"
        // );

        for (let i = 0; i < 1; i++) {
            for (let j = 0; j < 2; j++) {
                addColumnAgentGroup(
                    agentData,
                    1,
                    RADIUS * 1.5,
                    {
                        x: 10 - i * 6,
                        //x: 30,
                        z: 20 - j * 6,
                    },
                    {
                        x: 10 - i * 6,
                        //x: -20,
                        z: -999,
                    },
                    0.8,
                    "X"
                );

            }
        }

        // experiment border
        const boxGeometry1 = new THREE.BoxGeometry(31, 5, 1);
        const boxMaterial1 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
        const left = new THREE.Mesh(boxGeometry1, boxMaterial1);
        left.position.set(20, 2.5, 0);

        wallData.push({
            center: new THREE.Vector3(20, 2.5, 0),
            depth: 1, // along z-axis
            width: 31, // along x-axis
            base: new THREE.Vector3(20 - 31/2, 2.5, 0),
            tip: new THREE.Vector3(20 + 31/2, 2.5, 0),
        });

        const boxGeometry2 = new THREE.BoxGeometry(1, 5, 50);
        const boxMaterial2 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
        const bottom = new THREE.Mesh(boxGeometry2, boxMaterial2);
        //bottom.position.set(-15, 2.5, 0);
        bottom.position.set(-15, 2.5, -25);

        wallData.push({
            center: new THREE.Vector3(-15, 2.5, -25),
            depth: 50, // along z-axis
            width: 1, // along x-axis
            base: new THREE.Vector3(-15, 2.5, -25 - 50/2),
            tip: new THREE.Vector3(-15, 2.5, -25 + 50 / 2),
        });

        const boxGeometry3 = new THREE.BoxGeometry(1, 5, 50);
        const boxMaterial3 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
        const top = new THREE.Mesh(boxGeometry3, boxMaterial3);
        top.position.set(5, 2.5, -25);
        wallData.push({
            center: new THREE.Vector3(5, 2.5, -25),
            depth: 50, // along z-axis
            width: 1, // along x-axis
            base: new THREE.Vector3(5, 2.5, -25 - 50/2),
            tip: new THREE.Vector3(5, 2.5, -25 + 50/2),
        });

        const boxGeometry4 = new THREE.BoxGeometry(31, 5, 1);
        const boxMaterial4 = new THREE.MeshBasicMaterial({ color: 0x000f26 });
        const right = new THREE.Mesh(boxGeometry4, boxMaterial4);
        right.position.set(-30, 2.5, 0);
        wallData.push({
            center: new THREE.Vector3(-30, 2.5, 0),
            depth: 1, // along z-axis
            width: 31, // along x-axis
            base: new THREE.Vector3(-30 - 31 / 2, 2.5, 0),
            tip: new THREE.Vector3(-30 + 31 / 2, 2.5, 0),
        });

        scene.add(left);
        scene.add(bottom);
        scene.add(top);
        scene.add(right);

        parameters.wallData = wallData;

    }

//********************************************************* Trying scenario 3 *********************************************************

    function tryingScenario_Bilas_2p5() {

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 1.5,
            {
                x: 20,
                z: 2,
            },
            {
                x: -20,
                z: 2,
            },
            0.8,
            "X"
        );
        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 1.5,
            {
                x: 20,
                z: 2,
            },
            {
                x: -20,
                z: 2,
            },
            0.8,
            "X"
        );
        //-----------------------------------------------

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 4,
            {
                x: -20,
                z: 1
            },
            {
                x: 20,
                z: 1
            },
            0.8,
            "X"
        );

        addColumnAgentGroup(
            agentData,
            1,
            RADIUS * 4,
            {
                x: -20,
                z: 1
            },
            {
                x: 20,
                z: 1
            },
            0.8,
            "X"
        );



        // for(let i =0;i<agentData.length;i++){
        //
        //     if (i % 2 !== 0){
        //         agentData[i].agent.material = blueAgentMaterial;
        //
        //     }
        // }
    }



    function addColumnAgentGroup(
        agentData,
        numAgents,
        spacing,
        startPos,
        goalPos,
        velocityMagnitude,
        direction
    ) {
        let i = 0;
        let initalIdx = agentData.length;
        let dx = 0,
            dz = 0,
            vx = 0,
            vz = 0;
        let distanceToGoal = PHY.distance(
            startPos.x,
            startPos.z,
            goalPos.x,
            goalPos.z
        );
        vx = (velocityMagnitude * (goalPos.x - startPos.x)) / distanceToGoal;
        vz = (velocityMagnitude * (goalPos.z - startPos.z)) / distanceToGoal;

        if (direction == "X") {
            dx = spacing;
        } else if (direction == "Z") {
            dz = spacing;
        }
        while (i < numAgents) {
            agentData.push({
                index: i + initalIdx,
                x: startPos.x + dx * i,
                y: 2.0,
                z: startPos.z + dz * i,
                goal_x: goalPos.x + dx * i,
                goal_y: 0.0,
                goal_z: goalPos.z + dz * i,
                vx: vx,
                vy: 0.0,
                vz: vz,
                v_pref: Math.sqrt(vx * vx + vz * vz),
                radius: RADIUS,
                invmass: 0.5,
                colliding: false,
                group_id: 1,

                best: [],
                grad: {
                    'x':0,
                    'z':0,
                    'mx':1,
                    'mz':1,
                    's':0,
                    'dx':0,
                    'dz':0
                }
            });
            i += 1;
        }
    }

    let i = 0;
    let deltaSpacing = 3;
    let startX, startY, goalX, goalY;
    startX = -25;
    goalX = -25;
    startY = -20;
    goalY = 20;
    world.distanceConstraints = [];

    // defaultScenario();
    // testScenario()
    // testCrossScenario();58
    // testCrossWithDiagnoScenario();
    // testHallwayScenario();
    // testCrossScenario();
    // circleScenario();

    // tryingScenario_Bilas_1();
    //   tryingScenario_Bilas_1_debug();
    // tryingScenario_sphere_Bilas_1_debug();

    // tryingScenario_Bilas_2();
    // tryingScenario_Bilas_3();
      tryingScenario_Bilas_3_debug();
    // tryingScenario_Bilas_2p5();

    let agentGeom, agentMaterial, agent;
    let spotLight, spotLightTarget;
    let agentPointGeom, agentPointMaterial, agentPoint;

    agentData.forEach(function (item, index) {
        //agentGeom = new THREE.CylinderGeometry(item.radius, 1, 4, 16);
        // agentGeom = new THREE.CapsuleGeometry(item.radius, 2 * item.radius, 4, 8);
        agentGeom = new THREE.SphereGeometry( item.radius, 32, 16 );



        if (index % 2 !== 0){
            agentMaterial = new THREE.MeshLambertMaterial({
                color: 0x00ff00,
            });
        }else {
            agentMaterial = new THREE.MeshLambertMaterial({
                color: 0x0000ff,
            });
        }


        agent = new THREE.Mesh(agentGeom, agentMaterial);
        agent.castShadow = true;
        agent.receiveShadow = true;
        agent.userData = {
            index: item.index,
        };
        agent.rotateX(Math.PI / 2);
        // agent.rotateZ(Math.PI / 2);
        scene.add(agent);



        agentPointGeom = new THREE.CapsuleGeometry(item.radius, 2 * item.radius, 4, 8);
        agentPointMaterial = new THREE.MeshLambertMaterial({
            color: 0xffff00,
        });
        agentPoint = new THREE.Mesh(agentPointGeom, agentPointMaterial);
        agentPoint.castShadow = true;
        agentPoint.receiveShadow = true;
        // agent.rotateZ(Math.PI / 2);
        agentPoint.position.set(999, 2.5, 999);

        scene.add(agentPoint);

        // arrow for [x,z]
        let dir = new THREE.Vector3( 0, 1, 0 );
        let origin = agent.position;
        let length = 10;
        let hex = 0xffff00;

        let arrowHelper = new THREE.ArrowHelper( dir, origin, length, hex );
        // scene.add( arrowHelper );
        arrows.push(arrowHelper);

        // arrow for x and z
        let xdir = new THREE.Vector3( 0, 1, 0 );
        let xorigin = agent.position;
        let xhex = 0xffff00;
        let xarrowHelper = new THREE.ArrowHelper( xdir, xorigin, 1, xhex );
        scene.add( xarrowHelper );
        xarrows.push(xarrowHelper);

        let zdir = new THREE.Vector3( 0, 1, 0 );
        let zorigin = agent.position;
        let zhex = 0x00ffff;
        let zarrowHelper = new THREE.ArrowHelper( zdir, zorigin, 1, zhex );
        scene.add( zarrowHelper );
        zarrows.push(zarrowHelper);

        // velocity  indicator
        let g_dir = new THREE.Vector3( 1, 0, 0 );
        let g_origin = agent.position;
        let g_length = 5;
        let g_hex = 0x0000ff;

        let g_arrowHelper = new THREE.ArrowHelper( g_dir, g_origin, g_length, g_hex );
        // scene.add( g_arrowHelper );
        g_arrows.push(g_arrowHelper);




        // -----------------
        //adding spotlight code
        spotLight = new THREE.SpotLight(0xffffff);
        spotLight.position.set(item.x, item.y + 6, item.z);
        spotLight.shadow.mapSize.width = 1024;
        spotLight.shadow.mapSize.height = 1024;
        spotLight.shadow.camera.near = 500;
        spotLight.shadow.camera.far = 4000;
        spotLight.shadow.camera.fov = 30;
        spotLight.intensity = 0.4;
        spotLight.angle = Math.PI / 8;
        spotLightTarget = new THREE.Object3D();
        scene.add(spotLightTarget);
        spotLight.target = spotLightTarget;
        scene.add(spotLight);
        spotLights[item.index] = spotLight;
        // ----------------
        item.agent = agent;
        item.agentPoint = agentPoint;
        pickableObjects.push(agent);
        // pickableObjects.push(agent);

    });
    window.addEventListener("resize", onWindowResize);
    window.addEventListener("mousedown", mouseDown, false);
    // window.addEventListener("mousemove", mouseMove, false);

    // parameters setting
    for (let i =0; i<agentData.length;i++){
        parameters.best.push(
            Array(agentData.length).fill(null)
        )
    }

}

function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}

function mouseMove(event) {
    event.preventDefault();
    if (selected != null) {
        mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
        mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;
        raycaster.setFromCamera(mouse, camera);
        var intersects = raycaster.intersectObject(grid, false);
        for (let i = 0; i < intersects.length; i++) {

            agentData.forEach(function (member) {
                if (selected != null && member.index === selected) {
                    member.goal_x = intersects[i].point.x;
                    member.goal_z = intersects[i].point.z;
                    // ring.position.x = intersects[i].point.x;
                    // ring.position.z = intersects[i].point.z;
                }
            });
            break;
        }
    }
}

function mouseDown(event) {
    mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
    mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;
    raycaster.setFromCamera(mouse, camera);



    let intersects = raycaster.intersectObjects(pickableObjects, false);
    for (let i = 0; i < intersects.length; i++) {
        /* TODO finish this part as
         */
        let selectedUserData = intersects[i].object.userData;
        // console.log(agentData[selectedUserData.index]);
        console.log(agentData[selectedUserData.index].grad);
        console.log(selectedUserData.index);

        if (!selected.includes(selectedUserData.index)){
            selected.push(selectedUserData.index);
        }
        if (selected.length > 2){
            selected.shift();

            agentData.forEach(function (a){
                a.agentPoint.position.set(999, 2.5, 999);
            });
        }
        break;
    }
}

function render() {
    renderer.render(scene, camera);
}

function animate() {
    // console.log(mouse.x, mouse.y);
    requestAnimationFrame(animate);
    PHY.step(RADIUS, agentData, world, scene, parameters);


    // console.log(parameters.best);
    agentData.forEach(function (member, index) {
        // prevent agents from leaving the walls
        // if (member.x > 23) {
        //   member.x = 23;
        // }
        // if (member.x < -23) {
        //   member.x = -23;
        // }
        // if (member.z > 23) {
        //   member.z = 23;
        // }
        // if (member.z < -23) {
        //   member.z = -23;
        // }

        member.agent.position.x = member.x;
        member.agent.position.y = member.y;
        member.agent.position.z = member.z;

        // if(member.best !== null){
        //   member.agentPoint.position.x = member.best.x;
        //   member.agentPoint.position.y = member.y;
        //   member.agentPoint.position.z = member.best.z;
        // }


        const dx = member.goal_x - member.x;
        const dz = member.goal_z - member.z;
        member.agent.rotation.z = Math.atan2(dz, dx);

        member.agent.material = redAgentMaterial;

        if (member.colliding) {
            member.agent.material = greenAgentMaterial;
        }
        member.colliding = false;

        if (selected.length >0) {

            if (selected.length > 1){
                // visualize best

                let first = selected[0];
                let second = selected[1];

                let pair_best_points = parameters.best[first][second];

                let [first_best, second_best] = pair_best_points;
                agentData[first].agentPoint.position.x = first_best.x;
                agentData[first].agentPoint.position.z = first_best.z;

                agentData[second].agentPoint.position.x = second_best.x;
                agentData[second].agentPoint.position.z = second_best.z;
            }

            let include = selected.includes(member.index);
            if (include){
                member.agent.material = blueAgentMaterial;
            }


        }
        /* TODO finish this part as


             */
        spotLights[member.index].position.set(
            member.x - member.vx,
            member.y - member.vy,
            member.z - member.vz
        );
        spotLights[member.index].target.position.x = member.x;
        spotLights[member.index].target.position.y = member.y;
        spotLights[member.index].target.position.z = member.z;


        // visualizeXZMagnitude(member, index);
        // visualizeMagnitude(member, index)
        // visualizeVelocity(member, index);


    });
    renderer.render(scene, camera);
    stats.update();
}
init();
render();
animate();

// below are utiities
function visualizeXZMagnitude(member, index) {
    if (xarrows.length > 0) {

        xarrows[index].position.x = member.x;
        xarrows[index].position.y = member.y;
        xarrows[index].position.z = member.z;

        // console.log(member.grad);
        let direction = new THREE.Vector3(0, 1, 0);
        if (member.grad.x !== 0) {
            direction = new THREE.Vector3(member.grad.x, 0, 0);
        }

        xarrows[index].setDirection(direction.normalize());
        xarrows[index].setLength(direction.length() * 10);
    }

    if (zarrows.length > 0) {

        zarrows[index].position.x = member.x;
        zarrows[index].position.y = member.y;
        zarrows[index].position.z = member.z;

        // console.log(member.grad);
        let direction = new THREE.Vector3(0, 1, 0);
        if (member.grad.z !== 0) {
            direction = new THREE.Vector3(0, 0, member.grad.z);
        }

        zarrows[index].setDirection(direction.normalize());
        zarrows[index].setLength(direction.length() * 10);
    }

}
function visualizeMagnitude(member, index){
    if (arrows.length>0){

        arrows[index].position.x = member.x;
        arrows[index].position.y = member.y;
        arrows[index].position.z = member.z;

        // console.log(member.grad);
        let direction = new THREE.Vector3(0, 1, 0);
        if (member.grad.x !== 0 && member.grad.z !==0){
            direction = new THREE.Vector3(member.grad.x, 0, member.grad.z);
        }

        arrows[index].setDirection(direction.normalize());
        arrows[index].setLength(direction.length()*10);
    }
}


function visualizeVelocity(member, index){
    if (g_arrows.length>0){

        g_arrows[index].position.x = member.x;
        g_arrows[index].position.y = member.y;
        g_arrows[index].position.z = member.z;


        let direction = new THREE.Vector3(member.vx, 0, member.vz);


        g_arrows[index].setDirection(direction.normalize());
        g_arrows[index].setLength(direction.length()*10);
    }
}

