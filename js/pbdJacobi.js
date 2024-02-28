export function distance_2d(x1, y1, x2, y2) {
      return Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

export function distance(x1, y1, z1, x2, y2,z2) {
      return Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
}

//export function step(RADIUS,particles,world) {
export function step(particles, constraints) {
  /*  -----------------------  */
  /*  TODO modify lines below  */
  /*  -----------------------  */
  function distanceConstraint(agent_i,agent_j, desiredDistance,kScaler)
  {
    const agentCentroidDist = distance(agent_i.px, agent_i.py, agent_i.pz, 
                                       agent_j.px, agent_j.py, agent_j.pz);
    const agentDist = agentCentroidDist - desiredDistance;
    const dir_x = (agent_j.px- agent_i.px)/agentCentroidDist; 
    const dir_y = (agent_j.py- agent_i.py)/agentCentroidDist; 
    const dir_z = (agent_j.pz- agent_i.pz)/agentCentroidDist;
    const agent_i_scaler = kScaler*agent_i.invmass/(agent_i.invmass+agent_j.invmass) * agentDist;
    const agent_j_scaler = kScaler*agent_j.invmass/(agent_i.invmass+agent_j.invmass) * agentDist;
    if(Math.abs(agentDist) > epsilon )
    {
        agent_i.deltaX += agent_i_scaler * dir_x;
        agent_i.deltaY += agent_i_scaler * dir_y;
        agent_i.deltaZ += agent_i_scaler * dir_z;
        agent_i.deltaCtr +=1;
        agent_j.deltaX += -agent_j_scaler * dir_x;
        agent_j.deltaY += -agent_j_scaler * dir_y;
        agent_j.deltaZ += -agent_j_scaler * dir_z;
        agent_j.deltaCtr +=1;
    } 
  }


  function summarizeDeltas()
  {
    let item;
    for (let particleIndex in particles) {
      item = particles[particleIndex]; 
      if(item.deltaCtr>0)
      {
        item.px +=  item.deltaX/item.deltaCtr;
        item.py +=  item.deltaY/item.deltaCtr; 
        item.pz +=  item.deltaZ/item.deltaCtr; 
      }
    }
  }

  function resetDeltas()
  {
    let item;
    for (let particleIndex in particles) {
      item = particles[particleIndex];
      item.deltaX = 0;
      item.deltaY = 0;
      item.deltaZ = 0;
      item.deltaCtr = 0;
    }
  }

  function pointConstraint(agent_i,point)
  {
    const agentCentroidDist = distance(agent_i.px, agent_i.py, agent_i.pz, 
                                       point.x, point.y, point.z);
    const agentDist = agentCentroidDist;
    const dir_x = (point.x- agent_i.px)/agentCentroidDist; 
    const dir_y = (point.y- agent_i.py)/agentCentroidDist; 
    const dir_z = (point.z- agent_i.pz)/agentCentroidDist; 
    const agent_i_scaler = 1.0 * agentDist;
    if(Math.abs(agentDist) > epsilon )
    {
        agent_i.px +=  agent_i_scaler * dir_x
        agent_i.py +=  agent_i_scaler * dir_y
        agent_i.pz +=  agent_i_scaler * dir_z
    } 
  }
  const epsilon = 0.0001;
  const timestep = 0.01;
  const ITERNUM =14;
  for (let particleIndex in particles) {
    let item = particles[particleIndex]; 
  //particles.forEach(function (item) {
    //item.vy = -9.1*timestep + item.vy
    item.vx = 0.999*item.vx;
    item.vy = 0.999*item.vy;
    item.vz = 0.999*item.vz;
    item.px = item.x + timestep*item.vx; 
    item.pz = item.z + timestep*item.vz; 
    item.py = item.y + timestep*item.vy;
    item.deltaCtr = 0;
    item.deltaX = 0;
    item.deltaY = 0;
    item.deltaZ = 0;
  }

  let pbdIters = 0;
  const kStiffness = 0.99; // [0,1]
  var point, agent_a, agent_b, desDistance, i,j, idx = 0;
  /*
  while(pbdIters<ITERNUM)
  {
      resetDeltas();
      idx = 0;
      const kScaler = 1.0-(1.0-kStiffness)**(1.0/(pbdIters+1));
      while(idx < constraints.distanceConstraints.length)
      {
          desDistance = constraints.distanceConstraints[idx].distance;
          agent_a = particles[constraints.distanceConstraints[idx].idx_a]
          agent_b = particles[constraints.distanceConstraints[idx].idx_b]    
          distanceConstraint(agent_a,agent_b, desDistance,kScaler);
          idx+=1; 
      } 
      idx = 0;
      while(idx < constraints.pointConstraints.length)
      {
          agent_a = particles[constraints.pointConstraints[idx].idx_a]
          point = constraints.pointConstraints[idx].point;   
          pointConstraint(agent_a,point);
          idx+=1;
      
      }

      summarizeDeltas();

    pbdIters+=1;
  }
  */

  for (let particleIndex in particles) {
    let item = particles[particleIndex]; 
  //particles.forEach(function (item) {
    item.vx = (item.px-item.x)/timestep;
    item.vz = (item.pz-item.z)/timestep;
    item.vy = (item.py-item.y)/timestep;
    item.x = item.px; 
    item.z = item.pz; 
    item.y = item.py; 
  }

}
