import * as THREE from 'three';
import * as PHY from 'simplePhysics';
import {
    OrbitControls
} from "three/addons/controls/OrbitControls.js";


let renderer, scene, camera;
let world = {
    x: 80,
    z: 80
};
let agentData = [];
let pickableObjects = [];
let selected = null;
let mouse = new THREE.Vector2();
const raycaster = new THREE.Raycaster();
let grid,ring;

let spotLights = {};
let topTextures = {};
let topTexture;
const RADIUS = 1;
const blueAgentMaterial = new THREE.MeshLambertMaterial({
    color: 0x0000ff
});
const redAgentMaterial = new THREE.MeshLambertMaterial({
    color: 0xff0000
});


init();
render();


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
    camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 1, 1000);

    camera.position.set(-67.26, 54.07, -3.77);
    camera.rotation.order = 'YXZ';
    camera.rotation.y = -1.6267;
    camera.rotation.x = -0.46;

    // controls
    const controls = new OrbitControls(camera, renderer.domElement);
    controls.addEventListener('change', render);
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
    const texture = loader.load('resources/OIP.jpg');
    texture.wrapS = THREE.RepeatWrapping;
    texture.wrapT = THREE.RepeatWrapping;
    texture.magFilter = THREE.NearestFilter;
    const repeats = 40 / 32;
    texture.repeat.set(repeats, repeats);

    
    topTexture = loader.load('resources/triangle2.png');
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
    grid.rotation.order = 'YXZ';
    grid.rotation.y = -Math.PI / 2;
    grid.rotation.x = -Math.PI / 2;
    scene.add(grid);

    const ringGeometry = new THREE.RingGeometry( 1, 3, 12 );
    const ringMaterial = new THREE.MeshBasicMaterial( { color: 0xffff00, side: THREE.DoubleSide } );
    ring = new THREE.Mesh( ringGeometry, ringMaterial );
    scene.add( ring );
    ring.rotation.x = -Math.PI / 2;
    ring.position.y+=0.01;

    function addColumnAgentGroup(agentData, numAgents, spacing,
        startPos, goalPos,
        velocityMagnitude, direction) {
        let i = 0;
        let initalIdx = agentData.length;
        let dx = 0,
            dz = 0,
            vx = 0,
            vz = 0;
        let distanceToGoal = PHY.distance(startPos.x, startPos.z,
            goalPos.x, goalPos.z);
        vx = velocityMagnitude * (goalPos.x - startPos.x) / distanceToGoal;
        vz = velocityMagnitude * (goalPos.z - startPos.z) / distanceToGoal;

        if (direction == "X") {
            dx = spacing;
        } else if (direction == "Z") {
            dz = spacing;
        }
        while (i < numAgents) {
            agentData.push({
                index: i + initalIdx,
                x: startPos.x + dx * i,
                prev_xs: [],
                y: 2.0,
                z: startPos.z + dz * i,
                goal_x: goalPos.x + dx * i,
                goal_y: 0.0,
                goal_z: goalPos.z + dz * i,
                vx: vx,
                vy: 0.0,
                vz: vz,
                v_pref: Math.sqrt(vx*vx + vz*vz),
                radius: RADIUS,
                invmass: 0.5,
                group_id: 1
            });
            i += 1;
        }
    }


    let i = 0;
    let deltaSpacing = 3;
    let startX, startY, goalX, goalY;
    startX = -25;
    goalX = -25;
    startY = -20
    goalY = 20;
    world.distanceConstraints = [];

    addColumnAgentGroup(agentData, 3, RADIUS * 4, {
            x: 30,
            z: 25
        }, {
            x: -25,
            z: 25
        },
        0.6, "X", );

    addColumnAgentGroup(agentData, 4, RADIUS * 4, {
            x: 25,
            z: 20
        }, {
            x: -25,
            z: 20
        },
        0.7, "X", );

    addColumnAgentGroup(agentData, 4, RADIUS * 4, {
            x: 25,
            z: 10
        }, {
            x: -25,
            z: 10
        },
        0.8, "X", );

    addColumnAgentGroup(agentData, 4, RADIUS * 4, {
            x: 25,
            z: 6
        }, {
            x: -25,
            z: 6
        },
        0.8, "X", );


    let agnetGeometry, agentMaterial, agent;
    let spotLight, spotLightTarget;

    agentData.forEach(function(item) {
        agnetGeometry = new THREE.CylinderGeometry(item.radius, 1, 4, 16);
        agentMaterial = new THREE.MeshLambertMaterial({
            color: 0x00ff00
        }) ;
        //TODO add triangle on top of mesh and rotate 
        //see https://jsbin.com/zatejelika/1/edit?html,css,js,output

        agent = new THREE.Mesh(agnetGeometry, agentMaterial);
        agent.castShadow = true;
        agent.receiveShadow = true;
        agent.userData = {"index": item.index};
        scene.add(agent);
        // -----------------
        //adding spotlight code
        spotLight = new THREE.SpotLight( 0xffffff );
        spotLight.position.set( item.x, item.y+6, item.z);
        spotLight.shadow.mapSize.width = 1024;
        spotLight.shadow.mapSize.height = 1024;
        spotLight.shadow.camera.near = 500;
        spotLight.shadow.camera.far = 4000;
        spotLight.shadow.camera.fov = 30;
        spotLight.intensity=0.4;
        spotLight.angle = Math.PI/8;
        spotLightTarget = new THREE.Object3D();
        scene.add(spotLightTarget);
        spotLight.target = spotLightTarget;
        scene.add( spotLight );
        spotLights[item.index] = spotLight;
        // ----------------
        item.agent = agent;
        pickableObjects.push(agent);
    });
    window.addEventListener("resize", onWindowResize);
    window.addEventListener("mousedown", mouseDown, false);
    window.addEventListener("mousemove", mouseMove, false);
}    


function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(window.innerWidth, window.innerHeight);
}

function mouseMove(event) {
    event.preventDefault();
    if(selected!=null)
    {
        mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
        mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;
        raycaster.setFromCamera(mouse, camera);    
        var intersects =  raycaster.intersectObject(grid, false);
        for (var i = 0; i < intersects.length; i++) {
            agentData.forEach(function(member) {
                if(selected!=null&& member.index == selected)
                {
                    member.goal_x = intersects[i].point.x;
                    member.goal_z = intersects[i].point.z;
                    ring.position.x = intersects[i].point.x;
                    ring.position.z = intersects[i].point.z;
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
    selected=null;
    var intersects =  raycaster.intersectObjects(pickableObjects, false);
    for (var i = 0; i < intersects.length; i++) {
        selected = intersects[i].object.userData.index;
        break;
    }   
}

function render() {
    renderer.render(scene, camera);
}


function animate() {
    requestAnimationFrame(animate);
    PHY.step(RADIUS, agentData, world)
    agentData.forEach(function(member) {
        member.agent.position.x = member.x;
        member.agent.position.y = member.y;
        member.agent.position.z = member.z;
        //member.prev_xs.push({"x": member.x, "z": member.z});
        member.prev_xs.push(new THREE.Vector3(member.x, member.y+0.01,member.z));
        if(member.prev_xs.length > 10)
        {
            if(member.prev_xs.length > 350)
            {
                member.prev_xs.shift(); 
            }
            // Create a sine-like wave
            const splineCurve = new THREE.CatmullRomCurve3(member.prev_xs);
            const splinePoints = splineCurve.getPoints( 20 );
            const splinGeometry = new THREE.BufferGeometry().setFromPoints( splinePoints );
            const splineMaterial = new THREE.LineBasicMaterial( { color: 0xffffff } );
            // Create the final object to add to the scene
            const splineObject = new THREE.Line( splinGeometry, splineMaterial );
            scene.add(splineObject);
            if("line" in member)
            {
                scene.remove(member.line);
            }
            member.line = splineObject;
        }
        member.agent.material = redAgentMaterial;
        if(selected!=null&& member.index == selected)
        {
            member.agent.material = blueAgentMaterial;
        }
        //
        //if (member.group_id != 0) {
        //    member.agent.material = redAgentMaterial;
        //}
        //member.agent.material[1].map.rotation = Math.atan2(member.vy, member.vx);
        spotLights[member.index].position.set( member.x - member.vx, member.y - member.vy, member.z -member.vz);
        spotLights[member.index].target.position.x = member.x;// + member.vx; 
        spotLights[member.index].target.position.y = member.y;// + member.vy; 
        spotLights[member.index].target.position.z = member.z;// + member.vz; 

    
    });
    renderer.render(scene, camera);
};

animate();