% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PPR1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
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
%   see P3PPR1A0_convert_par2_MPV_fixb.m

% Output:
% tauX [3x1]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:37
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PPR1A0_invdyn_para_pf_mdp(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR1A0_invdyn_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PPR1A0_invdyn_para_pf_mdp: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PPR1A0_invdyn_para_pf_mdp: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR1A0_invdyn_para_pf_mdp: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PPR1A0_invdyn_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR1A0_invdyn_para_pf_mdp: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR1A0_invdyn_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR1A0_invdyn_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [6 1]), ...
  'P3PPR1A0_invdyn_para_pf_mdp: MDP has to be [6x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:37:32
% EndTime: 2019-05-03 14:37:32
% DurationCPUTime: 0.20s
% Computational Cost: add. (381->71), mult. (610->106), div. (0->0), fcn. (420->8), ass. (0->54)
t252 = xDDP(2) - g(2);
t253 = xDDP(1) - g(1);
t263 = xP(3);
t254 = sin(t263);
t255 = cos(t263);
t264 = koppelP(3,2);
t267 = koppelP(3,1);
t240 = t254 * t267 + t255 * t264;
t243 = -t254 * t264 + t255 * t267;
t256 = xDP(3) ^ 2;
t260 = xDDP(3);
t278 = t240 * t260 + t256 * t243 - t253;
t265 = koppelP(2,2);
t268 = koppelP(2,1);
t241 = t254 * t268 + t255 * t265;
t244 = -t254 * t265 + t255 * t268;
t277 = t241 * t260 + t256 * t244 - t253;
t266 = koppelP(1,2);
t269 = koppelP(1,1);
t242 = t254 * t269 + t255 * t266;
t245 = -t254 * t266 + t255 * t269;
t276 = t242 * t260 + t256 * t245 - t253;
t275 = -t256 * t240 + t243 * t260 + t252;
t274 = -t256 * t241 + t244 * t260 + t252;
t273 = -t256 * t242 + t245 * t260 + t252;
t257 = legFrame(3,3);
t246 = sin(t257);
t249 = cos(t257);
t218 = t278 * t246 + t275 * t249;
t258 = legFrame(2,3);
t247 = sin(t258);
t250 = cos(t258);
t219 = t277 * t247 + t274 * t250;
t259 = legFrame(1,3);
t248 = sin(t259);
t251 = cos(t259);
t220 = t276 * t248 + t273 * t251;
t234 = t246 * t267 - t249 * t264;
t235 = t246 * t264 + t249 * t267;
t236 = t247 * t268 - t250 * t265;
t237 = t247 * t265 + t250 * t268;
t238 = t248 * t269 - t251 * t266;
t239 = t248 * t266 + t251 * t269;
t272 = (t254 * t238 + t239 * t255) * t220 + (t254 * t236 + t237 * t255) * t219 + (t254 * t234 + t235 * t255) * t218;
t271 = t249 * t218 + t250 * t219 + t251 * t220;
t270 = -t246 * t218 - t247 * t219 - t248 * t220;
t233 = -t254 * t260 - t255 * t256;
t232 = -t254 * t256 + t255 * t260;
t231 = t254 * t252 + t255 * t253;
t230 = t255 * t252 - t254 * t253;
t223 = t273 * t248 - t276 * t251;
t222 = t274 * t247 - t277 * t250;
t221 = t275 * t246 - t278 * t249;
t1 = [t270 * MDP(1) + (t249 * t221 + t250 * t222 + t251 * t223 + t270) * MDP(2) + t233 * MDP(4) - t232 * MDP(5) + (-t254 * t230 + t255 * t231) * MDP(6); t271 * MDP(1) + (t246 * t221 + t247 * t222 + t248 * t223 + t271) * MDP(2) + t232 * MDP(4) + t233 * MDP(5) + (t255 * t230 + t254 * t231) * MDP(6); t272 * MDP(1) + ((t238 * t255 - t239 * t254) * t223 + (t236 * t255 - t237 * t254) * t222 + (t234 * t255 - t235 * t254) * t221 + t272) * MDP(2) + t260 * MDP(3) + t230 * MDP(4) - t231 * MDP(5);];
tauX  = t1;
