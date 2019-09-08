% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PPR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [6x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3PPR1G1P1A0_convert_par2_MPV_fixb.m

% Output:
% MMX [3x3]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:37
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMX = P3PPR1G1P1A0_inertia_para_pf_mdp(xP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR1G1P1A0_inertia_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR1G1P1A0_inertia_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR1G1P1A0_inertia_para_pf_mdp: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR1G1P1A0_inertia_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR1G1P1A0_inertia_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'P3PPR1G1P1A0_inertia_para_pf_mdp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:37:32
% EndTime: 2019-05-03 14:37:32
% DurationCPUTime: 0.11s
% Computational Cost: add. (159->58), mult. (285->92), div. (0->0), fcn. (266->8), ass. (0->43)
t322 = xP(3);
t317 = sin(t322);
t318 = cos(t322);
t319 = legFrame(3,3);
t314 = cos(t319);
t320 = legFrame(2,3);
t315 = cos(t320);
t321 = legFrame(1,3);
t316 = cos(t321);
t330 = t314 ^ 2 + t315 ^ 2 + t316 ^ 2;
t311 = sin(t319);
t312 = sin(t320);
t313 = sin(t321);
t331 = t311 ^ 2 + t312 ^ 2 + t313 ^ 2;
t334 = (t330 + t331) * MDP(2) + (t317 ^ 2 + t318 ^ 2) * MDP(6);
t323 = koppelP(3,2);
t326 = koppelP(3,1);
t292 = t311 * t326 - t314 * t323;
t293 = t311 * t323 + t314 * t326;
t285 = t317 * t292 + t293 * t318;
t324 = koppelP(2,2);
t327 = koppelP(2,1);
t294 = t312 * t327 - t315 * t324;
t295 = t312 * t324 + t315 * t327;
t287 = t317 * t294 + t295 * t318;
t325 = koppelP(1,2);
t328 = koppelP(1,1);
t296 = t313 * t328 - t316 * t325;
t297 = t313 * t325 + t316 * t328;
t289 = t317 * t296 + t297 * t318;
t298 = -t317 * t326 - t318 * t323;
t299 = -t317 * t327 - t318 * t324;
t300 = -t317 * t328 - t318 * t325;
t301 = -t317 * t323 + t318 * t326;
t302 = -t317 * t324 + t318 * t327;
t303 = -t317 * t325 + t318 * t328;
t333 = t289 * (-t300 * t313 + t303 * t316) + t287 * (-t299 * t312 + t302 * t315) + t285 * (-t298 * t311 + t301 * t314);
t332 = t285 * t314 + t287 * t315 + t289 * t316;
t329 = -t285 * t311 - t287 * t312 - t289 * t313;
t290 = t296 * t318 - t297 * t317;
t288 = t294 * t318 - t295 * t317;
t286 = t292 * t318 - t293 * t317;
t1 = [t331 * MDP(1) + t334; (-t314 * t311 - t315 * t312 - t316 * t313) * MDP(1); t330 * MDP(1) + t334; t329 * MDP(1) + (t286 * t314 + t288 * t315 + t290 * t316 + t329) * MDP(2) - t317 * MDP(4) - t318 * MDP(5); t332 * MDP(1) + (t286 * t311 + t288 * t312 + t290 * t313 + t332) * MDP(2) + t318 * MDP(4) - t317 * MDP(5); t333 * MDP(1) + (t290 * (t300 * t316 + t303 * t313) + t288 * (t299 * t315 + t302 * t312) + t286 * (t298 * t314 + t301 * t311) + t333) * MDP(2) + MDP(3);];
%% Postprocessing: Reshape Output
% From vec2symmat_3_matlab.m
res = [t1(1), t1(2), t1(4); t1(2), t1(3), t1(5); t1(4), t1(5), t1(6);];
MMX  = res;
