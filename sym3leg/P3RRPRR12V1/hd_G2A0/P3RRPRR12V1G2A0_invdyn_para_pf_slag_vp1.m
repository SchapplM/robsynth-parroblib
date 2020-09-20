% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR12V1G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:07
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:05:23
% EndTime: 2020-08-06 19:05:31
% DurationCPUTime: 7.15s
% Computational Cost: add. (31194->498), mult. (47145->818), div. (6867->6), fcn. (38043->18), ass. (0->322)
t208 = pkin(1) + pkin(2);
t197 = cos(qJ(2,3));
t205 = xDP(3);
t206 = xDP(2);
t207 = xDP(1);
t273 = t197 * t208;
t191 = sin(qJ(2,3));
t317 = qJ(3,3) * t191;
t149 = t273 + t317;
t192 = sin(qJ(1,3));
t198 = cos(qJ(1,3));
t233 = -pkin(4) * t192 + t149 * t198;
t143 = 0.1e1 / t149;
t210 = 0.1e1 / qJ(3,3);
t303 = t143 * t210;
t285 = t191 * t192;
t351 = pkin(4) * t198;
t140 = qJ(3,3) * t285 + t351;
t185 = legFrame(3,2);
t162 = sin(t185);
t165 = cos(t185);
t178 = t197 ^ 2;
t282 = t192 * t208;
t283 = t191 * t208;
t297 = t162 * qJ(3,3);
t85 = (t165 * t282 - t297) * t178 + (t140 * t165 + t162 * t283) * t197 + t297;
t291 = t165 * qJ(3,3);
t88 = (-t162 * t282 - t291) * t178 + (-t140 * t162 + t165 * t283) * t197 + t291;
t52 = (t197 * t205 * t233 + t206 * t88 + t207 * t85) * t303;
t334 = t208 * t52;
t199 = cos(qJ(2,2));
t271 = t199 * t208;
t193 = sin(qJ(2,2));
t319 = qJ(3,2) * t193;
t150 = t271 + t319;
t194 = sin(qJ(1,2));
t200 = cos(qJ(1,2));
t232 = -pkin(4) * t194 + t150 * t200;
t144 = 0.1e1 / t150;
t212 = 0.1e1 / qJ(3,2);
t302 = t144 * t212;
t281 = t193 * t194;
t350 = pkin(4) * t200;
t141 = qJ(3,2) * t281 + t350;
t186 = legFrame(2,2);
t163 = sin(t186);
t166 = cos(t186);
t179 = t199 ^ 2;
t278 = t194 * t208;
t279 = t193 * t208;
t295 = t163 * qJ(3,2);
t86 = (t166 * t278 - t295) * t179 + (t141 * t166 + t163 * t279) * t199 + t295;
t289 = t166 * qJ(3,2);
t89 = (-t163 * t278 - t289) * t179 + (-t141 * t163 + t166 * t279) * t199 + t289;
t53 = (t199 * t205 * t232 + t206 * t89 + t207 * t86) * t302;
t333 = t208 * t53;
t201 = cos(qJ(2,1));
t269 = t201 * t208;
t195 = sin(qJ(2,1));
t321 = qJ(3,1) * t195;
t151 = t269 + t321;
t196 = sin(qJ(1,1));
t202 = cos(qJ(1,1));
t231 = -pkin(4) * t196 + t151 * t202;
t145 = 0.1e1 / t151;
t214 = 0.1e1 / qJ(3,1);
t301 = t145 * t214;
t277 = t195 * t196;
t349 = pkin(4) * t202;
t142 = qJ(3,1) * t277 + t349;
t187 = legFrame(1,2);
t164 = sin(t187);
t167 = cos(t187);
t180 = t201 ^ 2;
t274 = t196 * t208;
t275 = t195 * t208;
t293 = t164 * qJ(3,1);
t87 = (t167 * t274 - t293) * t180 + (t142 * t167 + t164 * t275) * t201 + t293;
t287 = t167 * qJ(3,1);
t90 = (-t164 * t274 - t287) * t180 + (-t142 * t164 + t167 * t275) * t201 + t287;
t54 = (t201 * t205 * t231 + t206 * t90 + t207 * t87) * t301;
t332 = t208 * t54;
t209 = qJ(3,3) ^ 2;
t362 = 2 * rSges(3,3);
t365 = qJ(3,3) * t362 + t209;
t211 = qJ(3,2) ^ 2;
t364 = qJ(3,2) * t362 + t211;
t213 = qJ(3,1) ^ 2;
t363 = qJ(3,1) * t362 + t213;
t361 = m(1) * rSges(1,1);
t204 = m(2) * rSges(2,2);
t360 = m(2) * rSges(2,3);
t359 = (rSges(3,2) * m(3));
t82 = (-t192 * t205 + (-t162 * t206 + t165 * t207) * t198) * t143;
t358 = pkin(4) * t82;
t83 = (-t194 * t205 + (-t163 * t206 + t166 * t207) * t200) * t144;
t357 = pkin(4) * t83;
t84 = (-t196 * t205 + (-t164 * t206 + t167 * t207) * t202) * t145;
t356 = pkin(4) * t84;
t203 = pkin(1) + rSges(3,1);
t182 = -qJ(3,3) - rSges(3,3);
t355 = m(3) * t182;
t183 = -qJ(3,2) - rSges(3,3);
t354 = m(3) * t183;
t184 = -qJ(3,1) - rSges(3,3);
t353 = m(3) * t184;
t352 = m(3) * t203;
t348 = rSges(3,2) * t191;
t347 = rSges(3,2) * t193;
t346 = rSges(3,2) * t195;
t345 = t143 * t82;
t344 = t144 * t83;
t343 = t145 * t84;
t342 = t178 * t82;
t341 = t179 * t83;
t306 = t233 * t210;
t115 = t149 * t192 + t351;
t316 = qJ(3,3) * t197;
t224 = -t283 + t316;
t98 = -t115 * t162 - t165 * t224;
t328 = t210 * t98;
t97 = t115 * t165 - t162 * t224;
t329 = t210 * t97;
t67 = t205 * t306 + t206 * t328 + t207 * t329;
t340 = t182 * t67;
t305 = t232 * t212;
t116 = t150 * t194 + t350;
t318 = qJ(3,2) * t199;
t225 = -t279 + t318;
t100 = -t116 * t163 - t166 * t225;
t315 = t100 * t212;
t99 = t116 * t166 - t163 * t225;
t325 = t212 * t99;
t68 = t205 * t305 + t206 * t315 + t207 * t325;
t339 = t183 * t68;
t304 = t231 * t214;
t117 = t151 * t196 + t349;
t320 = qJ(3,1) * t201;
t226 = -t275 + t320;
t102 = -t117 * t164 - t167 * t226;
t313 = t102 * t214;
t101 = t117 * t167 - t164 * t226;
t314 = t101 * t214;
t69 = t205 * t304 + t206 * t313 + t207 * t314;
t338 = t184 * t69;
t22 = -t67 + t334;
t337 = t191 * t22;
t23 = -t68 + t333;
t336 = t193 * t23;
t24 = -t69 + t332;
t335 = t195 * t24;
t331 = t210 * t85;
t330 = t210 * t88;
t327 = t212 * t86;
t326 = t212 * t89;
t324 = t214 * t87;
t323 = t214 * t90;
t322 = t84 * t180;
t240 = -rSges(2,2) * t360 + Icges(2,6) - Icges(3,6);
t130 = -rSges(3,2) * t355 + t240;
t139 = rSges(2,1) * t360 + rSges(3,2) * t352 - Icges(3,4) - Icges(2,5);
t103 = t130 * t197 - t139 * t191;
t312 = t103 * t210;
t131 = -rSges(3,2) * t354 + t240;
t104 = t131 * t199 - t139 * t193;
t311 = t104 * t212;
t132 = -rSges(3,2) * t353 + t240;
t105 = t132 * t201 - t139 * t195;
t310 = t105 * t214;
t217 = rSges(2,2) ^ 2;
t219 = rSges(2,1) ^ 2;
t234 = (t217 + t219) * m(2) + Icges(3,2) + Icges(2,3);
t215 = rSges(3,3) ^ 2;
t222 = (pkin(1) ^ 2);
t235 = t215 + t222 + (2 * pkin(1) + rSges(3,1)) * rSges(3,1);
t106 = (t235 + t365) * m(3) + t234;
t309 = t106 * t210;
t107 = (t235 + t364) * m(3) + t234;
t308 = t107 * t212;
t108 = (t235 + t363) * m(3) + t234;
t307 = t108 * t214;
t160 = m(2) * rSges(2,1) + t352;
t300 = t160 * t197;
t299 = t160 * t199;
t298 = t160 * t201;
t296 = t162 * t198;
t294 = t163 * t200;
t292 = t164 * t202;
t290 = t165 * t198;
t288 = t166 * t200;
t286 = t167 * t202;
t284 = t191 * t197;
t280 = t193 * t199;
t276 = t195 * t201;
t272 = t197 * t210;
t270 = t199 * t212;
t268 = t201 * t214;
t267 = t203 * t210;
t266 = t203 * t212;
t265 = t203 * t214;
t264 = -Icges(2,1) - Icges(3,1);
t263 = rSges(3,2) ^ 2 + t215;
t262 = -2 * t359;
t261 = rSges(3,3) - t203;
t260 = rSges(3,3) + t203;
t259 = m(3) * t348;
t258 = m(3) * t347;
t257 = m(3) * t346;
t256 = m(3) * t340;
t255 = m(3) * t339;
t254 = m(3) * t338;
t253 = 0.2e1 * t284;
t252 = 0.2e1 * t280;
t251 = 0.2e1 * t276;
t250 = m(3) * t267;
t249 = m(3) * t266;
t248 = m(3) * t265;
t247 = t198 * t348;
t246 = t200 * t347;
t245 = t202 * t346;
t244 = t233 * t272;
t243 = t232 * t270;
t242 = t231 * t268;
t241 = -rSges(2,1) * t204 + Icges(2,4) - Icges(3,5);
t239 = -t222 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t154 = t204 + t355;
t155 = t204 + t354;
t156 = t204 + t353;
t238 = t210 * t259;
t237 = t212 * t258;
t236 = t214 * t257;
t136 = g(1) * t165 - g(2) * t162;
t230 = g(3) * t198 + t136 * t192;
t137 = g(1) * t166 - g(2) * t163;
t229 = g(3) * t200 + t137 * t194;
t138 = g(1) * t167 - g(2) * t164;
t228 = g(3) * t202 + t138 * t196;
t227 = Icges(2,2) + Icges(3,3) + (-t217 + t219) * m(2) + t264;
t223 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (rSges(2,3) ^ 2 + t217) * m(2) + Icges(1,3) - t264;
t220 = pkin(4) ^ 2;
t190 = xDDP(1);
t189 = xDDP(2);
t188 = xDDP(3);
t176 = g(3) * t361;
t153 = m(1) * rSges(1,2) - t359 - t360;
t152 = t153 * g(3);
t135 = g(1) * t164 + g(2) * t167;
t134 = g(1) * t163 + g(2) * t166;
t133 = g(1) * t162 + g(2) * t165;
t129 = -t184 * t352 + t241;
t128 = -t183 * t352 + t241;
t127 = -t182 * t352 + t241;
t114 = -(qJ(3,1) + t260) * (qJ(3,1) + t261) * m(3) + t227;
t113 = -(qJ(3,2) + t260) * (qJ(3,2) + t261) * m(3) + t227;
t112 = -(qJ(3,3) + t260) * (qJ(3,3) + t261) * m(3) + t227;
t81 = t84 ^ 2;
t80 = t83 ^ 2;
t79 = t82 ^ 2;
t78 = t195 * t356;
t77 = t193 * t357;
t76 = t191 * t358;
t75 = t114 * t180 + t129 * t251 + (t263 + t363) * m(3) + t223;
t74 = t113 * t179 + t128 * t252 + (t263 + t364) * m(3) + t223;
t73 = t112 * t178 + t127 * t253 + (t263 + t365) * m(3) + t223;
t72 = (t304 + (-rSges(3,2) * t277 - t203 * t242) * t145) * m(3);
t71 = (t305 + (-rSges(3,2) * t281 - t203 * t243) * t144) * m(3);
t70 = (t306 + (-rSges(3,2) * t285 - t203 * t244) * t143) * m(3);
t66 = -t231 * t248 + (-t105 * t196 + t108 * t242) * t145;
t65 = -t232 * t249 + (-t104 * t194 + t107 * t243) * t144;
t64 = -t233 * t250 + (-t103 * t192 + t106 * t244) * t143;
t63 = (t314 + (t167 * t245 - t87 * t265) * t145) * m(3);
t62 = (t325 + (t166 * t246 - t86 * t266) * t144) * m(3);
t61 = (t329 + (t165 * t247 - t85 * t267) * t143) * m(3);
t60 = (t313 + (-t164 * t245 - t90 * t265) * t145) * m(3);
t59 = (t315 + (-t163 * t246 - t89 * t266) * t144) * m(3);
t58 = (t328 + (-t162 * t247 - t88 * t267) * t143) * m(3);
t57 = t231 * t236 + (t105 * t242 - t196 * t75) * t145;
t56 = t232 * t237 + (t104 * t243 - t194 * t74) * t144;
t55 = t233 * t238 + (t103 * t244 - t192 * t73) * t143;
t51 = t54 ^ 2;
t50 = t53 ^ 2;
t49 = t52 ^ 2;
t42 = -t101 * t248 + (t105 * t286 + t87 * t307) * t145;
t41 = -t99 * t249 + (t104 * t288 + t86 * t308) * t144;
t40 = -t97 * t250 + (t103 * t290 + t85 * t309) * t143;
t39 = -t102 * t248 + (-t105 * t292 + t90 * t307) * t145;
t38 = -t100 * t249 + (-t104 * t294 + t89 * t308) * t144;
t37 = -t98 * t250 + (-t103 * t296 + t88 * t309) * t143;
t36 = t129 * t54;
t35 = t128 * t53;
t34 = t127 * t52;
t33 = t101 * t236 + (t75 * t286 + t87 * t310) * t145;
t32 = t99 * t237 + (t74 * t288 + t86 * t311) * t144;
t31 = t97 * t238 + (t73 * t290 + t85 * t312) * t143;
t30 = t102 * t236 + (-t75 * t292 + t90 * t310) * t145;
t29 = t100 * t237 + (-t74 * t294 + t89 * t311) * t144;
t28 = t98 * t238 + (-t73 * t296 + t88 * t312) * t143;
t27 = t78 + t332;
t26 = t77 + t333;
t25 = t76 + t334;
t21 = (-t54 * t320 + t335) * pkin(4) + ((qJ(3,1) + t208) * (-qJ(3,1) + t208) * t180 + qJ(3,1) * t208 * t251 + t213 + t220) * t84;
t20 = (-t53 * t318 + t336) * pkin(4) + ((qJ(3,2) + t208) * (-qJ(3,2) + t208) * t179 + qJ(3,2) * t208 * t252 + t211 + t220) * t83;
t19 = (-t52 * t316 + t337) * pkin(4) + ((qJ(3,3) + t208) * (-qJ(3,3) + t208) * t178 + t208 * qJ(3,3) * t253 + t209 + t220) * t82;
t18 = (-t356 + (t226 + t320) * t54 + (0.2e1 * t69 - t332) * t195) * t343;
t17 = (-t357 + (t225 + t318) * t53 + (0.2e1 * t68 - t333) * t193) * t344;
t16 = (-t358 + (t224 + t316) * t52 + (0.2e1 * t67 - t334) * t191) * t345;
t15 = (((-t213 + t239) * t54 + t208 * t69) * t54 + t27 * t69 + (t226 * t54 * pkin(4) - t21) * t84) * t214;
t14 = (((-t211 + t239) * t53 + t208 * t68) * t53 + t26 * t68 + (t225 * t53 * pkin(4) - t20) * t83) * t212;
t13 = (((-t209 + t239) * t52 + t208 * t67) * t52 + t25 * t67 + (t224 * t52 * pkin(4) - t19) * t82) * t210;
t12 = -t21 * t268 * t343 + ((-(t78 + t24) * t269 + (pkin(4) * t322 - t335) * qJ(3,1)) * t54 + (t201 * t27 + t54 * t321) * t69) * t301;
t11 = -t20 * t270 * t344 + ((-(t77 + t23) * t271 + (pkin(4) * t341 - t336) * qJ(3,2)) * t53 + (t199 * t26 + t53 * t319) * t68) * t302;
t10 = -t19 * t272 * t345 + ((-(t76 + t22) * t273 + (pkin(4) * t342 - t337) * qJ(3,3)) * t52 + (t197 * t25 + t52 * t317) * t67) * t303;
t9 = t12 * t352 - t18 * t257 + (-t15 - t81 * t203 * t276 - (t180 * t81 - t51 - t81) * t184 + t135 * t201 - t228 * t195) * m(3);
t8 = t11 * t352 - t17 * t258 + (-t14 - t80 * t203 * t280 - (t179 * t80 - t50 - t80) * t183 + t134 * t199 - t229 * t193) * m(3);
t7 = t10 * t352 - t16 * t259 + (-t13 - t79 * t203 * t284 - (t178 * t79 - t49 - t79) * t182 + t133 * t197 - t230 * t191) * m(3);
t6 = -t105 * t18 - t108 * t12 + (-t135 * t160 + t228 * t156) * t201 + (t135 * t156 + t228 * t160) * t195 + (t15 * t203 - 0.2e1 * t54 * t338) * m(3) + (t114 * t276 - 0.2e1 * t129 * t180 + t129) * t81;
t5 = -t104 * t17 - t107 * t11 + (-t134 * t160 + t229 * t155) * t199 + (t134 * t155 + t229 * t160) * t193 + (t14 * t203 - 0.2e1 * t53 * t339) * m(3) + (t113 * t280 - 0.2e1 * t128 * t179 + t128) * t80;
t4 = -t103 * t16 - t106 * t10 + (-t133 * t160 + t230 * t154) * t197 + (t133 * t154 + t230 * t160) * t191 + (t13 * t203 - 0.2e1 * t52 * t340) * m(3) + (t112 * t284 - 0.2e1 * t127 * t178 + t127) * t79;
t3 = -t75 * t18 - t105 * t12 - 0.4e1 * (-t36 - t254 / 0.2e1) * t322 - (t139 * t54 + t69 * t262) * t54 * t201 + 0.2e1 * t84 * (-t36 - t254) + t152 * t202 + (g(3) * t298 + t176) * t196 + (-t15 * t359 - 0.2e1 * t84 * (t114 * t54 - t69 * t352) * t201 - t132 * t51 - g(3) * t156 * t196) * t195 + ((t156 * t195 - t298 - t361) * t202 + t153 * t196) * t138;
t2 = -t74 * t17 - t104 * t11 - 0.4e1 * (-t35 - t255 / 0.2e1) * t341 - (t139 * t53 + t68 * t262) * t53 * t199 + 0.2e1 * t83 * (-t35 - t255) + t152 * t200 + (g(3) * t299 + t176) * t194 + (-t14 * t359 - 0.2e1 * t83 * (t113 * t53 - t68 * t352) * t199 - t131 * t50 - g(3) * t155 * t194) * t193 + ((t155 * t193 - t299 - t361) * t200 + t153 * t194) * t137;
t1 = -t73 * t16 - t103 * t10 - 0.4e1 * (-t34 - t256 / 0.2e1) * t342 - (t139 * t52 + t67 * t262) * t52 * t197 + 0.2e1 * t82 * (-t34 - t256) + t152 * t198 + (g(3) * t300 + t176) * t192 + (-t13 * t359 - 0.2e1 * t82 * (t112 * t52 - t67 * t352) * t197 - t130 * t49 - g(3) * t154 * t192) * t191 + ((t154 * t191 - t300 - t361) * t198 + t153 * t192) * t136;
t43 = [t9 * t314 + t7 * t329 + t8 * t325 - m(4) * g(1) + (t63 * t313 + t62 * t315 + t61 * t328) * t189 + (t63 * t304 + t62 * t305 + t61 * t306) * t188 + (t63 * t314 + t62 * t325 + t61 * t329 + m(4)) * t190 + ((t33 * t286 + t42 * t324) * t190 + (-t33 * t292 + t42 * t323) * t189 + (-t196 * t33 + t42 * t242) * t188 + t3 * t286 + t6 * t324) * t145 + ((t32 * t288 + t41 * t327) * t190 + (-t32 * t294 + t41 * t326) * t189 + (-t194 * t32 + t41 * t243) * t188 + t2 * t288 + t5 * t327) * t144 + ((t31 * t290 + t40 * t331) * t190 + (-t31 * t296 + t40 * t330) * t189 + (-t192 * t31 + t40 * t244) * t188 + t1 * t290 + t4 * t331) * t143; t8 * t315 + t9 * t313 + t7 * t328 - m(4) * g(2) + (t60 * t314 + t59 * t325 + t58 * t329) * t190 + (t60 * t304 + t59 * t305 + t58 * t306) * t188 + (t60 * t313 + t59 * t315 + t58 * t328 + m(4)) * t189 + ((t30 * t286 + t39 * t324) * t190 + (-t30 * t292 + t39 * t323) * t189 + (-t196 * t30 + t39 * t242) * t188 - t3 * t292 + t6 * t323) * t145 + ((t29 * t288 + t38 * t327) * t190 + (-t29 * t294 + t38 * t326) * t189 + (-t194 * t29 + t38 * t243) * t188 - t2 * t294 + t5 * t326) * t144 + ((t28 * t290 + t37 * t331) * t190 + (-t28 * t296 + t37 * t330) * t189 + (-t192 * t28 + t37 * t244) * t188 - t1 * t296 + t4 * t330) * t143; t7 * t306 + t8 * t305 + t9 * t304 - m(4) * g(3) + (t72 * t314 + t71 * t325 + t70 * t329) * t190 + (t72 * t313 + t71 * t315 + t70 * t328) * t189 + (t72 * t304 + t71 * t305 + t70 * t306 + m(4)) * t188 + ((t57 * t286 + t66 * t324) * t190 + (-t57 * t292 + t66 * t323) * t189 + (-t196 * t57 + t66 * t242) * t188 - t196 * t3 + t6 * t242) * t145 + ((t56 * t288 + t65 * t327) * t190 + (-t56 * t294 + t65 * t326) * t189 + (-t194 * t56 + t65 * t243) * t188 - t194 * t2 + t5 * t243) * t144 + ((t55 * t290 + t64 * t331) * t190 + (-t55 * t296 + t64 * t330) * t189 + (-t192 * t55 + t64 * t244) * t188 - t192 * t1 + t4 * t244) * t143;];
tauX  = t43;
