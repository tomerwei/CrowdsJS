import * as THREE from 'three';
import * as PHY from 'simplePhysics';
import {
    OrbitControls
} from "three/addons/controls/OrbitControls.js";

import Stats from 'three/addons/libs/stats.module.js';

let renderer, scene, camera;
let world = {
    x: 80,
    z: 80
};
let agentData = [];

// added variables
let agentGoalFlags = [];
let agentStartFlags = []; //0 is starting at the top, 1 is starting at the side
//

let pickableObjects = [];
let wallsData = [];
let selected = null;
let mouse = new THREE.Vector2();
const raycaster = new THREE.Raycaster();
let grid,ring;

//added variables
let goalArea = 6; //"radius" by which the goals are extended
let curgoalPos = {x: 0, z: 0};
let goalPos1 = {x: 0, z: 0}; //first goal, center of the intersection
let goalPos2 = {x: -50,z: 0}; //second goal, the "exit"
//

let spotLights = {};

let topTexture;
const RADIUS = 1;
const blueAgentMaterial = new THREE.MeshLambertMaterial({
    color: 0x0000ff
});
const redAgentMaterial = new THREE.MeshLambertMaterial({
    color: 0xff0000
});

    const stats = new Stats();
    document.body.appendChild(stats.dom)   

init();
render();


function init() {
	curgoalPos = goalPos1;
	
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

    function addColumnAgentGroup(agentData, numAgents, spacing, startPos, velocityMagnitude, direction) {
		console.log("b")
        let initalIdx = agentData.length;
        let vx = 0,
            vz = 0;
		let distanceToGoal = PHY.distance(startPos.x, startPos.z, curgoalPos.x, curgoalPos.z);
		vx = velocityMagnitude * (curgoalPos.x - startPos.x) / distanceToGoal;
		vz = velocityMagnitude * (curgoalPos.z - startPos.z) / distanceToGoal;

        for (let i = 0; i < numAgents; i++) {
			console.log("a")
            agentData.push({
                index: i + initalIdx,
                x: startPos.x,
                y: 2.0,
                z: startPos.z,
				goal_x: curgoalPos.x,
				goal_y: 0.0,
				goal_z: curgoalPos.z,
                vx: vx,
                vy: 0.0,
                vz: vz,
                v_pref: Math.sqrt(vx*vx + vz*vz),
                radius: RADIUS,
                invmass: 0.5
            });
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
	let num_of_agents_up = 100 //represents the agents that will start at the top hallway
	let num_of_agents_side = 100 //represents the agents that will start at the side hallway
	let up_x_max = 35;
	let up_x_min = 10;
	let up_z_max = 5;
	let up_z_min = -5;
	let side_x_max = 5;
	let side_x_min = -5;
	let side_z_max = 35;
	let side_z_min = 10;
	for (let i = 0; i < num_of_agents_up; i++) {
		addColumnAgentGroup(agentData, 1, RADIUS * 4, { //start position
            x: getRandomInt(up_x_min, up_x_max),
            z: getRandomInt(up_z_min, up_z_max)
        }, 10, "X", );
		agentGoalFlags.push(false);
	}
	for (let i = 0; i < num_of_agents_side; i++) {
		addColumnAgentGroup(agentData, 1, RADIUS * 4, { //start position
            x: getRandomInt(side_x_min, side_x_max),
            z: -getRandomInt(side_z_min, side_z_max)
        }, 10, "X", );
		agentGoalFlags.push(false);
	}
    let agnetGeometry, agentMaterial, agent;
    let spotLight, spotLightTarget;

    agentData.forEach(function(item) {
        agnetGeometry = new THREE.CylinderGeometry(item.radius, 1, 4, 16);
        agentMaterial = new THREE.MeshLambertMaterial({
            color: 0x00ff00
        }) ;

        agent = new THREE.Mesh(agnetGeometry, agentMaterial);
        agent.castShadow = true;
        agent.receiveShadow = true;
        agent.userData = {"index": item.index};
        scene.add(agent);

        item.agent = agent;
        pickableObjects.push(agent);
    });

//new walls created o make the shape of the hallway/intersection
    wallsData.push({
		//these represent position
                "x": -22.5,
                "y": 0,
                "z": -22.5,
				//these represent size of the wall
                "dx":30.0,
                "dy":15.0,
                "dz":30.0,
            });
    wallsData.push({
		//these represent position
                "x": 22.5,
                "y": 0,
                "z": -22.5,
				//these represent size of the wall
                "dx":30.0,
                "dy":15.0,
                "dz":30.0,
            });
	wallsData.push({
		//these represent position
                "x": 0,
                "y": 0,
                "z": 17.5,
				//these represent size of the wall
                "dx":75.0,
                "dy":15.0,
                "dz":22.5,
            });
			
    let wallGeometry, wall, wallMaterial;
    wallsData.forEach(function (item) {
        wallGeometry = new THREE.BoxGeometry(item.dx, item.dy, item.dz);
        wallMaterial = new THREE.MeshLambertMaterial({
            color: 0x00ff00
        }) ;

        wall = new THREE.Mesh(wallGeometry, wallMaterial);
        wall.castShadow = true;
        wall.receiveShadow = true;
        wall.userData = {"index": item.index};
        scene.add(wall);
        wall.position.x = item.x; 
        wall.position.y = item.y; 
        wall.position.z = item.z; 
        //pickableObjects.push(wall);
    });
	world['wallData'] = wallsData
	//world - data place that is everything that is NOT an agent
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
            /*            
            agentData.forEach(function(member) {
                if(selected!=null )
                {

                }
            });
            */
            break;
        }   
    }
}
function getRandomInt(min_int, max_int) {
	return Math.random() * (max_int - min_int) + min_int;
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
	let count = 0;
    agentData.forEach(function(member) {
        member.agent.position.x = member.x;
        member.agent.position.y = member.y;
        member.agent.position.z = member.z;
		//console.log(member.agent.position.x);
		//console.log(member.agent.position.z);
		if (member.agent.position.x <= 40 && member.agent.position.x >= -10 && member.agent.position.z <= goalPos1.z + goalArea && member.agent.position.z >= goalPos1.z - goalArea) {
			//console.log("first goal reached!")
			agentGoalFlags[count] = true;
			console.log(agentGoalFlags[count])
		}
        member.agent.material = redAgentMaterial;
        if(selected!=null&& member.index == selected)
        {
            member.agent.material = blueAgentMaterial;
        }
		if (member.agent.position.x <= -35) {
			agentGoalFlags[count] = false;
			let rand = getRandomInt(0, 100); //random number to determine which end of the intersection the agent/member will be teleported too.
			if (rand >= 50) {
				member.x = 39.5;
				member.z = getRandomInt(-5, 5);
				agentStartFlags[count] = 1;
			}
			else {
				member.x = getRandomInt(-5, 5);
				member.z = -39.5;
				agentStartFlags[count] = 0;
			}
		}
		count++;

    });
	for (let i = 0; i < agentGoalFlags.length; i++) {
		//console.log(agentGoalFlags[i])
		if (agentGoalFlags[i] == true) {
			agentData[i].goal_x = goalPos2.x;
			agentData[i].goal_z = agentData[i].z;
			//console.log("next goal time");
		}
		else {
			agentData[i].goal_x = goalPos1.x;
			agentData[i].goal_z = goalPos1.z;
			if (agentStartFlags[i] == 1) {
				agentData[i].goal_z = agentData[i].z;
			}
			if (agentStartFlags[i] == 0) {
				agentData[i].goal_x = agentData[i].x;
			}
		}
	}

    renderer.render(scene, camera);
    stats.update()
};

animate();