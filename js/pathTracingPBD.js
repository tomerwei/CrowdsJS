export function distance(x1, y1, x2, y2) {
      return Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

export function step(RADIUS,sceneEntities,world) {


  /*  -----------------------  */
  /*  TODO modify lines below  */
  /*  -----------------------  */
  function distanceConstraint(agent_i,agent_j, desiredDistance)
  {
    const agentCentroidDist = distance(agent_i.px, agent_i.pz, 
                agent_j.px, agent_j.pz );
    const agentDist = agentCentroidDist - desiredDistance;
    const dir_x = (agent_j.px- agent_i.px)/agentCentroidDist; 
    const dir_z = (agent_j.pz- agent_i.pz)/agentCentroidDist; 
    const agent_i_scaler = 0.1*agent_i.invmass/(agent_i.invmass+agent_j.invmass) * agentDist;
    const agent_j_scaler = 0.1*agent_j.invmass/(agent_i.invmass+agent_j.invmass) * agentDist;
    if(Math.abs(agentDist) > epsilon )
    {
        agent_i.px +=  agent_i_scaler * dir_x
        agent_i.pz +=  agent_i_scaler * dir_z
        agent_j.px += - agent_j_scaler * dir_x
        agent_j.pz += - agent_j_scaler * dir_z
    } 
  }

  function collisionConstraint(agent_i,agent_j)
  {
    const agentCentroidDist = distance(agent_i.px, agent_i.pz, 
                agent_j.px, agent_j.pz );
    const agentDist = agentCentroidDist - AGENTSIZE;
    const dir_x = (agent_j.px- agent_i.px)/agentCentroidDist; 
    const dir_z = (agent_j.pz- agent_i.pz)/agentCentroidDist;
    const agent_i_scaler = agent_i.invmass/(agent_i.invmass+agent_j.invmass) * agentDist
    const agent_j_scaler = agent_j.invmass/(agent_i.invmass+agent_j.invmass) * agentDist 
    if(agentDist < 0)
    {
        agent_i.px += agent_i_scaler * dir_x
        agent_i.pz += agent_i_scaler * dir_z
        agent_j.px += - agent_j_scaler * dir_x
        agent_j.pz += - agent_j_scaler * dir_z
    } 
  }


  function longRangeConstraint(agent_i,agent_j)
  {
    const dist = distance(agent_i.px, agent_i.pz, 
                agent_j.px, agent_j.pz );
    const radius_init = AGENTSIZE;
    const radius_sq = radius_init*radius_init;
    if (dist < radius_init)
    {
      radius_sq = (radius_init - dist) * (radius_init - dist)
    }
    const v_x = (1.0/time_delta) * ((agent_i.px - agent_i.x) - (agent_j.px - agent_j.x))
    const v_z = (1.0/time_delta) * ((agent_i.pz - agent_i.z) - (agent_j.pz - agent_j.z))
    const x0_x = agent_i.x - agent_j.x;
    const x0_z = agent_i.z - agent_j.z;


    const dir_x = (agent_j.px- agent_i.px)/dist; 
    const dir_z = (agent_j.pz- agent_i.pz)/dist;
    const agent_i_scaler = agent_i.invmass/(agent_i.invmass+agent_j.invmass) * agentDist
    const agent_j_scaler = agent_j.invmass/(agent_i.invmass+agent_j.invmass) * agentDist 
    if(agentDist < 0)
    {
        agent_i.px += agent_i_scaler * dir_x
        agent_i.pz += agent_i_scaler * dir_z
        agent_j.px += - agent_j_scaler * dir_x
        agent_j.pz += - agent_j_scaler * dir_z
    } 
/*
for p_i in positions:
        for j in range(particle_num_neighbors[p_i]):
            p_j = particle_neighbors[p_i, j]
            if p_j >= 0 and (group_number[p_j] != group_number[p_i] or not is_shape_matching_agent[p_i]):
                p_j = particle_neighbors[p_i, j]
                radius_init = agent_radii[p_j] + agent_radii[p_i]
                radius_sq_init = radius_init * radius_init
                pos_ji = positions[p_i] - positions[p_j]
                dist = pos_ji.norm()
                radius_sq = radius_sq_init
                if (dist < radius_init):
                    radius_sq = (radius_init - dist) * (radius_init - dist)
                v_ = (1.0/time_delta) * ((positions[p_i] - old_positions[p_i]) - (
                    positions[p_j] - old_positions[p_j]))
                x0_ = old_positions[p_i] - old_positions[p_j]
                v_sq = v_.norm_sqr()
                x_sqaure_vec = x0_ * x0_
                x_sq = x_sqaure_vec.x + x_sqaure_vec.y
                a = v_sq
                b = -1 * v_.dot(x0_)
                b_sq = b * b
                c = x_sq - radius_sq
                d_sq = b_sq - a * c
                d = ti.sqrt(d_sq)
                tao = (b - d) / a
                if (d_sq > 0.0 and ti.abs(a) > 0.0001 and tao > 0 and tao < C_TAU_MAX):
                    dv_i = 1.
                    dv_j = -1.
                    clamp_tao = ti.exp(-tao * tao / tao0)
                    c_tao = abs(tao - tao0)
                    v_x = v_[0]
                    v_y = v_[1]
                    x0 = x0_[0]
                    y0 = x0_[1]
                    x0_sq = x_sqaure_vec.x
                    y0_sq = x_sqaure_vec.y
                    x0 = x0_.x
                    y0 = x0_.y
                    grad_x_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_x * tao) -
                                                          (x0 + (v_y * x0 * y0 + v_x * (radius_sq - y0_sq)) / d)))
                    grad_y_i = 2 * c_tao * ((dv_i / a) * ((-2. * v_y * tao) -
                                                          (y0 + (v_x * x0 * y0 + v_y * (radius_sq - x0_sq)) / d)))
                    grad_x_j = -grad_x_i
                    grad_y_j = -grad_y_i
                    stiff = ti.exp(-tao * tao / tao0)
                    tao_sq = c_tao * c_tao
                    stiff = ti.exp(-tao * tao / 25.0)
                    combined_radii = agent_radii[p_i] + agent_radii[p_j]
                    inv_mass_i = agent_inv_mass[p_i]
                    inv_mass_j = agent_inv_mass[p_j]
                    w_i_coef = inv_mass_i / (inv_mass_i + inv_mass_j)
                    w_j_coef = inv_mass_j / (inv_mass_i + inv_mass_j)
                    s = stiff * tao_sq / (inv_mass_i *
                                          (grad_y_i * grad_y_i + grad_x_i * grad_x_i) +
                                          inv_mass_j *
                                          (grad_y_j * grad_y_j + grad_x_j * grad_x_j))
                    cur_deltas = ti.Vector([0.0, 0.0])
                    cur_deltas.x = s * w_i_coef * grad_x_i
                    cur_deltas.y = s * w_i_coef * grad_y_i
                    #cur_deltas = clamp(cur_deltas, C_MAX_ACCELERATION)
                    lengthV = cur_deltas.norm()
                    if lengthV > C_LONG_RANGE_CONSTRAINT:
                        mult = C_LONG_RANGE_CONSTRAINT / lengthV
                        cur_deltas = mult * cur_deltas
                    position_deltas[p_i].x += cur_deltas.x
                    position_deltas[p_i].y += cur_deltas.y
                    position_deltas_ctr[p_i] += 1






*/
  }


  function agentVelocityPlanner()  {
    sceneEntities.forEach(function (agent_i) {
      const distToGoal = distance(agent_i.x, agent_i.z, 
                        agent_i.goal_x, agent_i.goal_z );
      if(distToGoal > RADIUS)
      {
        const dir_x = (agent_i.goal_x- agent_i.x)/distToGoal; 
        const dir_z = (agent_i.goal_z- agent_i.z)/distToGoal;
        agent_i.vx = agent_i.v_pref * dir_x;
        agent_i.vz = agent_i.v_pref * dir_z;
      }
      agent_i.vx = 0.9999*agent_i.vx;
      agent_i.vz = 0.9999*agent_i.vz;      
    });
  }

  /*  -----------------------  */



  const AGENTSIZE = RADIUS * 2; 
  const epsilon = 0.0001;
  const timestep = 0.03;
  const ITERNUM =3;

  agentVelocityPlanner();

  sceneEntities.forEach(function (item) {
    item.px = item.x + timestep*item.vx; 
    item.pz = item.z + timestep*item.vz; 
    item.py = item.y + timestep*item.vy;
  });


  let pbdIters = 0;
  var agent_a, agent_b, desDistance, i,j, idx = 0;

  while(pbdIters<ITERNUM)
  {
      idx = 0;
      while(idx < world.distanceConstraints.length)
      {
          desDistance = world.distanceConstraints[idx].distance;
          agent_a = sceneEntities[world.distanceConstraints[idx].idx_a]
          agent_b = sceneEntities[world.distanceConstraints[idx].idx_b]         
          distanceConstraint(agent_a,agent_b, desDistance);
          idx+=1; 
      }  
      i=0;
      while(i<sceneEntities.length)
      {
          j=i+1;
          while(j<sceneEntities.length)
          {
            collisionConstraint(sceneEntities[i],sceneEntities[j])
            j+=1;
          }
          i+=1
      }
    pbdIters+=1;
  }


  sceneEntities.forEach(function (item) {
    item.vx = (item.px-item.x)/timestep;
    item.vz = (item.pz-item.z)/timestep;
    item.vy = (item.py-item.y)/timestep;
    item.x = item.px; 
    item.z = item.pz; 
    item.y = item.py; 

    if(item.x < -world.x/2)
    {
      item.x = -world.x/2;
    }
    else if(item.x > world.x/2)
    {
      item.x = world.x/2; 
    }
    if(item.z < -world.z/2)
    {
      item.z = -world.z/2;
    }
    else if(item.z > world.z/2)
    {
      item.z= world.z/2; 
    }
  });

}
