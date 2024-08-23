import * as THREE from "three";

export function distance(x1, y1, x2, y2) {
  return Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

// Function to calculate the angle between two vectors
function angleBetweenVectors_2(v1, v2) {
  // Calculate dot product
  const dotProduct = v1.x * v2.x + v1.z * v2.z;

  // Calculate magnitudes of vectors
  const magnitudeV1 = Math.sqrt(v1.x * v1.x + v1.z * v1.z);
  const magnitudeV2 = Math.sqrt(v2.x * v2.x + v2.z * v2.z);

  // Calculate angle in radians using dot product and magnitudes
  const angleRadians = Math.acos(dotProduct / ((magnitudeV1 * magnitudeV2)+0.0000001) );

  // Convert angle from radians to degrees
  const angleDegrees = angleRadians * 180 / Math.PI;

  return angleDegrees;
}
function getRandomNum(min_num, max_num) {
	return Math.random() * (max_num - min_num) + min_num;
}


export function step(RADIUS, sceneEntities, world, scene, customParams = {}) {

	let timestep = 0
	timestep = 0.03;

	const ITERNUM = 1; // 3
	const agentLength = RADIUS;
	const WallsData = world['wallsData'];

	// collision functions
	function rotateLineSegment(x1, y1, x2, y2, r) {
		// Calculate the center of the line segment
		const centerX = (x1 + x2) / 2;
		const centerY = (y1 + y2) / 2;

		// Translate the line segment so that its center is at the origin
		const x1p = x1 - centerX;
		const y1p = y1 - centerY;
		const x2p = x2 - centerX;
		const y2p = y2 - centerY;

		// Rotate the line segment about the origin
		const cosR = Math.cos(r);
		const sinR = Math.sin(r);
		const x1r = x1p * cosR - y1p * sinR;
		const y1r = x1p * sinR + y1p * cosR;
		const x2r = x2p * cosR - y2p * sinR;
		const y2r = x2p * sinR + y2p * cosR;

		// Translate the line segment back to its original position
		const newX1 = x1r + centerX;
		const newY1 = y1r + centerY;
		const newX2 = x2r + centerX;
		const newY2 = y2r + centerY;

		// Return the new endpoints of the line segment
		return [newX1, newY1, newX2, newY2];
	 }

	function getCapsuleBodyNormal(agent, agentLength, RADIUS, current_rotation) {
		let iCoords = rotateLineSegment(
			agent.x,
			agent.z + agentLength + RADIUS,
			agent.x,
			agent.z - agentLength - RADIUS,
			current_rotation
		);

		if(customParams.scenario==='bottleneck') {
			iCoords = rotateLineSegment(
				agent.z,
				agent.x + agentLength + RADIUS,
				agent.z,
				agent.x - agentLength - RADIUS,
				current_rotation
			);
		};
	  
		const aa = {
			tip: new THREE.Vector3(iCoords[0], 0, iCoords[1]),
			base: new THREE.Vector3(iCoords[2], 0, iCoords[3]),
		};
      
		// Calculate the slope of the line
		const dx = aa.base.x - aa.tip.x;
		const dz = aa.base.z - aa.tip.z  ;
		
		let leng = Math.sqrt(dx*dx + dz*dz);
			
		// Generate a normal vector (-dy, dx), perpendicular to the line
		let nx = -dz / leng;  
		let nz = dx / leng;
	  
		let normal_to_capsule = new THREE.Vector3(nx, 0, nz);
		agent.normal_to_capsule = normal_to_capsule;

		  
		let angle_capsule_normal_and_vel = angleBetweenVectors_2(agent.normal_to_capsule, agent.normal_to_capsule_prev);

		if(angle_capsule_normal_and_vel > 20){    // forcing normal vector to capsule body to be in the capsule's facing direction. 
			agent.normal_to_capsule = agent.normal_to_capsule_prev;
		}
		agent.normal_to_capsule_prev = agent.normal_to_capsule;
		  
		return agent.normal_to_capsule;
	}

	  /*  -----------------------  */
	  /*  TODO modify lines below  */
	  /*  -----------------------  */

	function collisionConstraint(agent_i,agent_j){
		const agentCentroidDist = distance(agent_i.px, agent_i.pz, agent_j.px, agent_j.pz );
		const agentDist_x = agentCentroidDist - ((agentLength + RADIUS) * 2);
		const agentDist_z = agentCentroidDist - (RADIUS * 2);
		const dir_x = (agent_j.px- agent_i.px)/agentCentroidDist; 
		const dir_z = (agent_j.pz- agent_i.pz)/agentCentroidDist;
		const agent_i_scaler = agent_i.invmass/(agent_i.invmass+agent_j.invmass) * agentDist_x
		const agent_j_scaler = agent_j.invmass/(agent_i.invmass+agent_j.invmass) * agentDist_x
		if(agentDist_x < 0) //since the agents aren't shaped normally, check for x and z distance seperately
		{
			console.log("COLLIDING");
			agent_i.px += agent_i_scaler * dir_x		
			agent_j.px += - agent_j_scaler * dir_x
		} 
		if (agentDist_z < 0) {
			console.log("COLLIDING");
			agent_j.pz += - agent_j_scaler * dir_z
			agent_i.pz += agent_i_scaler * dir_z
		}
	}

	function resolveObstacleCollision(agent, obstacle) {
		let obstacleX_TopLeft = obstacle.x - obstacle.dx/2;
		let obstacleX_BottomRight = obstacle.x + obstacle.dx/2;
		let obstacleZ_TopLeft = obstacle.z - obstacle.dz/2;
		let obstacleZ_BottomRight = obstacle.z + obstacle.dz/2;
		let normalDirX, normalDirZ;
		if (agent.px + RADIUS > obstacleX_TopLeft && agent.pz + RADIUS > obstacleZ_TopLeft && agent.px - RADIUS < obstacleX_BottomRight && agent.pz - RADIUS < obstacleZ_BottomRight) { //this if statement is never becoming true
			console.log("yes");
			normalDirX = (agent.px - (obstacleX_TopLeft + obstacleX_BottomRight)/2);
			normalDirZ = (agent.pz - (obstacleZ_TopLeft + obstacleZ_BottomRight)/2);
			let xDist = Math.abs(normalDirX);
			let zDist = Math.abs(normalDirZ);
			if (xDist > zDist) {
				normalDirZ /= zDist;
				normalDirX = 0.0;
			}
			else {
				normalDirZ = 0.0;
				normalDirX /= xDist;
			}
			let w_A = 1.0;
			let w_B = 0.0;
			agent.px += w_A / (w_A + w_B) * xDist * normalDirX;
			agent.pz += w_A / (w_A + w_B) * zDist * normalDirZ;
		}
	}
		
	function agentVelocityPlanner() {
		sceneEntities.forEach(function (agent_i) {
			const distToGoal = distance(
				agent_i.x,
				agent_i.z,
				agent_i.goal_x,
				agent_i.goal_z
			);
			if (distToGoal > RADIUS) {
				const dir_x = (agent_i.goal_x - agent_i.x) / distToGoal;
				const dir_z = (agent_i.goal_z - agent_i.z) / distToGoal;
				agent_i.vx = agent_i.v_pref * dir_x;
				agent_i.vz = agent_i.v_pref * dir_z;
			}
			agent_i.vx = 0.9999 * agent_i.vx;
			agent_i.vz = 0.9999 * agent_i.vz;
			if(customParams.scenario == 'dense_torso_like' )
			{
			  // ---------- for Torso Dense Crowd -------  START -----------------------------------------------
				// if(distToGoal < 5 * RADIUS)
				if(distToGoal < 3 * RADIUS && agent_i.index != 0)
				{
				  agent_i.vx = 0.01 * agent_i.vx;
				  agent_i.vz = 0.01 * agent_i.vz;        
				}
				//----------- for Torso Dense Crowd ------- END ----------------------------------------------
			}

		 
		});
	}

	function toVisualize_Capsule_normal(capsule_entity){
		let current_rotation = capsule_entity.agent.rotation.z;
		getCapsuleBodyNormal(capsule_entity, agentLength, RADIUS, current_rotation);
	}


	/*  -----------------------  */
	agentVelocityPlanner();
	
	sceneEntities.forEach(function (item) {
		item.px = item.x + timestep * item.vx;
		item.pz = item.z + timestep * item.vz;
		item.py = item.y + timestep * item.vy;
	});

	let pbdIters = 0;
	let isColliding;
	var agent_a,
		agent_b,
		desDistance,
		i,
		j,
		k,
		idx = 0;
	while (pbdIters < ITERNUM) {

		// clean previous accumulated gradient
		i = 0;
		while (i < sceneEntities.length) {
			j = i + 1;
			while (j < sceneEntities.length) {
				collisionConstraint(sceneEntities[i],sceneEntities[j]);
				sceneEntities[i].grad.x = 0;
				sceneEntities[i].grad.z = 0;
				sceneEntities[j].grad.x = 0;
				sceneEntities[j].grad.z = 0;

				sceneEntities[i].grad.dx = 0;
				sceneEntities[i].grad.dz = 0;
				sceneEntities[j].grad.dx = 0;
				sceneEntities[j].grad.dz = 0;
				j += 1;
			}
			k = 0;
			/*while (k<WallsData.length) {
				resolveObstacleCollision(sceneEntities[i], WallsData[k])
				k+=1;
			}*/
			console.log(WallsData);
			i += 1;
		}

		function rotationConstraint_V2(capsule_entity)
		{
			let capsuleBodyNormalVec = getCapsuleBodyNormal(capsule_entity, agentLength, RADIUS, capsule_entity.agent.rotation.z);
			let VelocityVec = new THREE.Vector3(capsule_entity.vx, 0 , capsule_entity.vz);
			// let capsuleCurToGoalVec = new THREE.Vector3(capsule_entity.goal_z - capsule_entity.z, 0, capsule_entity.goal_x - capsule_entity.x);
			let angleBodyNormalToGoalVec = angleBetweenVectors_2(capsuleBodyNormalVec, VelocityVec);
			let cur_orientation = capsule_entity.agent.rotation.z;

			//smooth the rotation speed while changing orientation
			if( (Math.abs(capsule_entity.agent.rotation.z - capsule_entity.nextOrientationInRadians) >= 0.08)  && (angleBodyNormalToGoalVec < 90)   )   //20 for rectangle
			{
				if(customParams.scenario == 'swap_Scenario') {
					capsule_entity.agent.rotation.z = capsule_entity.agent.rotation.z + capsule_entity.nextOrientationInRadians/200;        // 150  200
				}
			}

			// Rotate the velocity vector by 90 degrees in the 2D plane to get the perpendicular vector
			const perpendicularVector = new THREE.Vector3(-capsule_entity.vz, 0, capsule_entity.vx);
			// Compute the dot product
			const dotProduct = perpendicularVector.dot(capsuleBodyNormalVec);
			// Determine if the point is to the left or right of the velocity vector.
			const direction = dotProduct > 0 ? 'right' : 'left';

			let  next_orientation = 0;	
			let angleBodyNormalToGoalVecInRad = angleBodyNormalToGoalVec * (Math.PI / 180);

			// find the shortest-path rotation. 
			if(direction == 'right')
			{
				next_orientation = cur_orientation - angleBodyNormalToGoalVecInRad;
			}else{
				next_orientation = cur_orientation + angleBodyNormalToGoalVecInRad;
			}
		  

			if( customParams.orientation == 'front' && angleBodyNormalToGoalVec > 1 )
			{
				if( cur_orientation >= next_orientation )
				{
					capsule_entity.agent.rotation.z = cur_orientation - 0.2; //THIS CONTROLS THE rotation
				}else{
					capsule_entity.agent.rotation.z = cur_orientation + 0.2;
				}   
			}
		}
		i = 0;
		while (i < sceneEntities.length) {
			rotationConstraint_V2(sceneEntities[i]);
			i += 1;
		}
		i = 0;
		while (i < sceneEntities.length) {
			toVisualize_Capsule_normal(sceneEntities[i])
			i += 1;
		}
		pbdIters += 1;
	}
 

	sceneEntities.forEach(function (item) {
		item.vx = (item.px - item.x) / timestep;
		item.vz = (item.pz - item.z) / timestep;
		item.vy = (item.py - item.y) / timestep;

		item.x = item.px;
		item.z = item.pz;
		item.y = item.py;

// ------------------------- Start------------------ For swap_Scenario  -------------------------------------------------------------------------------------------------------------

	if( (customParams.scenario === 'swap_Scenario')  )
	{
		let dist3 = distance(item.x, item.z, item.goal_x, item.goal_z); //first goal
		if(dist3 < 10.0 )
		{
			item.goal_x = item.x;
			item.goal_z = -40;
		}

		let dist = distance(item.x, item.z, item.x, -40);  //second goal 
		if( dist < 5.0 ){
			let randomNum = getRandomNum(-10, 10);
			//console.log(item.s_Flag);
			if (item.s_Flag == 0) { 
				//console.log("yes");
				item.x = -35;
				item.z = randomNum;
				item.goal_x = 0;
				item.goal_z = randomNum;
			}
			else {
				//console.log("no");
				item.x = randomNum;
				item.z = 35;
				item.goal_x = randomNum;
				item.goal_z = 0;
			}
		}
	}

	if( (distance(item.x, item.z, item.goal_x, item.goal_z) < 1 ) && ( item.z > 15  || item.z < -15)  )
	{
		item.vx = 0;
		item.vy = 0;
		item.vz = 0;

		item.x = item.goal_x;
		item.z = item.goal_z;
		item.y = 0;

		if( item.index == 0 )
		{
			item.agent.rotation.z = 0; 
		}    
	}
// ------------------------- End------------------ For swap_Scenario  -------------------------------------------

	});

}
