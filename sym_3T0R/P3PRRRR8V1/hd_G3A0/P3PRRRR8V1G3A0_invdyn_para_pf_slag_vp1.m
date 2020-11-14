% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRRR8V1G3A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 17:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:15:48
% EndTime: 2020-08-06 17:15:55
% DurationCPUTime: 6.75s
% Computational Cost: add. (18822->396), mult. (58296->786), div. (7479->8), fcn. (61890->28), ass. (0->322)
t181 = cos(qJ(3,2));
t182 = cos(qJ(2,2));
t176 = sin(qJ(2,2));
t275 = t176 * t181;
t127 = pkin(2) * t275 - pkin(5) * t182;
t164 = sin(pkin(3));
t175 = sin(qJ(3,2));
t166 = cos(pkin(3));
t351 = pkin(2) * t166;
t109 = t127 * t164 + t175 * t351;
t370 = 0.1e1 / t109;
t321 = t370 / t181;
t179 = cos(qJ(3,3));
t180 = cos(qJ(2,3));
t174 = sin(qJ(2,3));
t277 = t174 * t179;
t126 = pkin(2) * t277 - pkin(5) * t180;
t173 = sin(qJ(3,3));
t108 = t126 * t164 + t173 * t351;
t371 = 0.1e1 / t108;
t322 = t371 / t179;
t169 = legFrame(1,2);
t152 = sin(t169);
t155 = cos(t169);
t122 = g(1) * t152 + g(2) * t155;
t125 = g(1) * t155 - g(2) * t152;
t165 = cos(pkin(6));
t149 = g(3) * t165;
t163 = sin(pkin(6));
t240 = -t125 * t163 - t149;
t381 = t122 * t164 + t240 * t166;
t168 = legFrame(2,2);
t151 = sin(t168);
t154 = cos(t168);
t121 = g(1) * t151 + g(2) * t154;
t124 = g(1) * t154 - g(2) * t151;
t242 = -t124 * t163 - t149;
t380 = t121 * t164 + t242 * t166;
t167 = legFrame(3,2);
t150 = sin(t167);
t153 = cos(t167);
t120 = g(1) * t150 + g(2) * t153;
t123 = g(1) * t153 - g(2) * t150;
t244 = -t123 * t163 - t149;
t379 = t120 * t164 + t244 * t166;
t177 = sin(qJ(3,1));
t183 = cos(qJ(3,1));
t375 = -rSges(3,1) * t183 + rSges(3,2) * t177;
t374 = -rSges(3,1) * t181 + rSges(3,2) * t175;
t373 = -rSges(3,1) * t179 + rSges(3,2) * t173;
t188 = xDP(3);
t189 = xDP(2);
t190 = xDP(1);
t235 = t150 * t189 - t153 * t190;
t286 = t166 * t174;
t211 = t164 * t179 + t173 * t286;
t278 = t173 * t180;
t88 = t211 * t163 - t165 * t278;
t91 = t163 * t278 + t211 * t165;
t64 = (t188 * t88 + t235 * t91) * t322;
t61 = t64 ^ 2;
t234 = t151 * t189 - t154 * t190;
t285 = t166 * t176;
t210 = t164 * t181 + t175 * t285;
t276 = t175 * t182;
t89 = t210 * t163 - t165 * t276;
t92 = t163 * t276 + t210 * t165;
t65 = (t188 * t89 + t234 * t92) * t321;
t62 = t65 ^ 2;
t233 = t152 * t189 - t155 * t190;
t162 = 0.1e1 / t183;
t178 = sin(qJ(2,1));
t273 = t178 * t183;
t248 = t164 * t273;
t284 = t166 * t177;
t184 = cos(qJ(2,1));
t287 = t164 * t184;
t372 = 0.1e1 / (-pkin(5) * t287 + (t248 + t284) * pkin(2));
t320 = t372 * t162;
t283 = t166 * t178;
t209 = t164 * t183 + t177 * t283;
t274 = t177 * t184;
t90 = t209 * t163 - t165 * t274;
t93 = t163 * t274 + t209 * t165;
t66 = (t188 * t90 + t233 * t93) * t320;
t63 = t66 ^ 2;
t147 = m(3) * rSges(3,1) * rSges(3,2) - Icges(3,4);
t369 = 0.2e1 * t147;
t187 = m(2) * rSges(2,1);
t368 = m(3) * rSges(3,3);
t338 = rSges(3,2) * t179;
t132 = rSges(3,1) * t173 + t338;
t94 = -t164 * t174 * t132 - t166 * t373;
t367 = m(3) * t94;
t337 = rSges(3,2) * t181;
t133 = rSges(3,1) * t175 + t337;
t95 = -t164 * t176 * t133 - t166 * t374;
t366 = m(3) * t95;
t336 = rSges(3,2) * t183;
t134 = rSges(3,1) * t177 + t336;
t96 = -t164 * t178 * t134 - t166 * t375;
t365 = m(3) * t96;
t198 = 0.1e1 / pkin(2);
t282 = t166 * t180;
t350 = pkin(2) * t179;
t82 = (-t163 * t174 + t165 * t282) * t350 + pkin(5) * (t163 * t180 + t165 * t286);
t85 = (t163 * t282 + t165 * t174) * t350 + (t163 * t286 - t165 * t180) * pkin(5);
t58 = (t188 * t85 + t235 * t82) * t198 * t322;
t364 = pkin(2) * t58;
t281 = t166 * t182;
t349 = pkin(2) * t181;
t83 = (-t163 * t176 + t165 * t281) * t349 + pkin(5) * (t163 * t182 + t165 * t285);
t86 = (t163 * t281 + t165 * t176) * t349 + (t163 * t285 - t165 * t182) * pkin(5);
t59 = (t188 * t86 + t234 * t83) * t198 * t321;
t363 = pkin(2) * t59;
t280 = t166 * t184;
t348 = pkin(2) * t183;
t84 = (-t163 * t178 + t165 * t280) * t348 + pkin(5) * (t163 * t184 + t165 * t283);
t87 = (t163 * t280 + t165 * t178) * t348 + (t163 * t283 - t165 * t184) * pkin(5);
t60 = (t188 * t87 + t233 * t84) * t198 * t320;
t362 = pkin(2) * t60;
t361 = pkin(5) * t64;
t360 = pkin(5) * t65;
t359 = pkin(5) * t66;
t194 = rSges(3,2) ^ 2;
t195 = rSges(3,1) ^ 2;
t136 = (-t194 + t195) * m(3) - Icges(3,1) + Icges(3,2);
t358 = t136 / 0.2e1;
t145 = rSges(3,2) * t368 - Icges(3,6);
t357 = -t145 / 0.4e1;
t146 = rSges(3,1) * t368 - Icges(3,5);
t356 = t146 / 0.4e1;
t355 = m(3) * t198;
t157 = t179 ^ 2;
t354 = pkin(2) * t157;
t159 = t181 ^ 2;
t353 = pkin(2) * t159;
t161 = t183 ^ 2;
t352 = pkin(2) * t161;
t347 = g(3) * t163;
t196 = pkin(5) ^ 2;
t197 = pkin(2) ^ 2;
t331 = t173 * t58;
t265 = pkin(2) * t331;
t346 = (-pkin(5) * t265 + (t157 * t197 + t196) * t64) * t64;
t330 = t175 * t59;
t264 = pkin(2) * t330;
t345 = (-pkin(5) * t264 + (t159 * t197 + t196) * t65) * t65;
t335 = t372 * t66;
t117 = -m(3) * t373 + t187;
t144 = m(2) * rSges(2,2) - t368;
t97 = t117 * t180 - t144 * t174;
t334 = t164 * t97;
t118 = -m(3) * t374 + t187;
t98 = t118 * t182 - t144 * t176;
t333 = t164 * t98;
t119 = -m(3) * t375 + t187;
t99 = t119 * t184 - t144 * t178;
t332 = t164 * t99;
t329 = t177 * t60;
t328 = t198 * t82;
t327 = t198 * t83;
t326 = t198 * t84;
t325 = t198 * t85;
t324 = t198 * t86;
t323 = t198 * t87;
t156 = m(1) + m(2) + m(3);
t319 = t371 * t156;
t318 = t370 * t156;
t128 = pkin(2) * t273 - pkin(5) * t184;
t110 = pkin(2) * t284 + t128 * t164;
t107 = 0.1e1 / t110;
t317 = t107 * t156;
t114 = -t145 * t179 - t146 * t173;
t316 = t114 * t198;
t115 = -t145 * t181 - t146 * t175;
t315 = t115 * t198;
t116 = -t145 * t183 - t146 * t177;
t314 = t116 * t198;
t312 = t120 * t166;
t310 = t121 * t166;
t308 = t122 * t166;
t304 = t136 * t173;
t303 = t136 * t175;
t302 = t136 * t177;
t272 = t194 + t195;
t141 = t272 * m(3) + Icges(3,3);
t301 = t141 * t198;
t300 = rSges(3,2) * t164 * t149;
t171 = xDDP(2);
t299 = t150 * t171;
t298 = t151 * t171;
t297 = t152 * t171;
t172 = xDDP(1);
t296 = t153 * t172;
t295 = t154 * t172;
t294 = t155 * t172;
t293 = t163 * t164;
t292 = t164 * t173;
t291 = t164 * t175;
t290 = t164 * t177;
t289 = t164 * t180;
t288 = t164 * t182;
t279 = t166 * t198;
t271 = t371 * t367;
t270 = t370 * t366;
t269 = t107 * t365;
t268 = t94 * t355;
t267 = t95 * t355;
t266 = t96 * t355;
t263 = pkin(2) * t329;
t262 = t173 * t361;
t261 = t175 * t360;
t260 = t177 * t359;
t259 = t371 * t334;
t258 = t370 * t333;
t257 = t107 * t332;
t256 = t150 * t322;
t255 = t153 * t322;
t254 = t151 * t321;
t253 = t154 * t321;
t252 = t152 * t320;
t251 = t155 * t320;
t250 = t164 * t277;
t249 = t164 * t275;
t243 = -t123 * t165 + t347;
t241 = -t124 * t165 + t347;
t239 = -t125 * t165 + t347;
t22 = t262 - t364;
t10 = (((t166 * t58 + t64 * t289) * t354 - (t265 - t361) * t250 + t166 * t22) * t64 - (-t58 * t289 + (-t157 * t166 + t173 * t250 + t166) * t64) * t364) * t322;
t13 = (t279 * t346 + (-t58 * t126 * t292 + t166 * (t58 * t354 - t262)) * t58) * t322;
t16 = (t22 * t364 - t346) * t371;
t191 = 0.2e1 * qJ(3,3);
t202 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,1) / 0.2e1 + Icges(3,2) / 0.2e1 + (0.2e1 * rSges(3,3) ^ 2 + t272) * m(3) / 0.2e1;
t73 = cos(t191) * t358 - t147 * sin(t191) + t202;
t4 = -t16 * t334 - t73 * t10 - t114 * t13 - 0.4e1 * ((t64 * t304 / 0.2e1 + t58 * t356) * t179 + t331 * t357 + (t157 - 0.1e1 / 0.2e1) * t64 * t147) * t58 + (-t117 * t379 - t243 * t144) * t180 - (t243 * t117 - t144 * t379) * t174;
t201 = t379 * t174 - t243 * t180;
t7 = -t16 * t367 - t114 * t10 - t141 * t13 + t61 * (t157 * t369 + t179 * t304 - t147) + (((t244 * t164 - t312) * rSges(3,1) + t201 * rSges(3,2)) * t179 + t173 * (t300 + (t123 * t293 + t312) * rSges(3,2) + t201 * rSges(3,1))) * m(3);
t232 = -t7 * t328 - t91 * t4;
t23 = t261 - t363;
t11 = (((t166 * t59 + t65 * t288) * t353 - (t264 - t360) * t249 + t166 * t23) * t65 + (t59 * t288 + (t159 * t166 - t175 * t249 - t166) * t65) * t363) * t321;
t14 = (t279 * t345 + (-t59 * t127 * t291 + t166 * (t59 * t353 - t261)) * t59) * t321;
t17 = (t23 * t363 - t345) * t370;
t192 = 0.2e1 * qJ(3,2);
t74 = cos(t192) * t358 - t147 * sin(t192) + t202;
t5 = -t17 * t333 - t74 * t11 - t115 * t14 - 0.4e1 * ((t65 * t303 / 0.2e1 + t59 * t356) * t181 + t330 * t357 + (t159 - 0.1e1 / 0.2e1) * t65 * t147) * t59 + (-t118 * t380 - t241 * t144) * t182 - (t241 * t118 - t144 * t380) * t176;
t200 = t380 * t176 - t241 * t182;
t8 = -t17 * t366 - t115 * t11 - t141 * t14 + t62 * (t159 * t369 + t181 * t303 - t147) + (((t242 * t164 - t310) * rSges(3,1) + t200 * rSges(3,2)) * t181 + t175 * (t300 + (t124 * t293 + t310) * rSges(3,2) + t200 * rSges(3,1))) * m(3);
t231 = -t8 * t327 - t92 * t5;
t24 = t260 - t362;
t12 = (((t166 * t60 + t66 * t287) * t352 - (t263 - t359) * t248 + t166 * t24) * t335 + (t60 * t287 + (t161 * t166 - t177 * t248 - t166) * t66) * t372 * t362) * t162;
t21 = -pkin(5) * t263 + (t161 * t197 + t196) * t66;
t15 = (t21 * t279 * t335 + (-t60 * t128 * t290 + t166 * (t60 * t352 - t260)) * t107 * t60) * t162;
t18 = (t21 * t66 - t24 * t362) * t372;
t193 = 0.2e1 * qJ(3,1);
t75 = cos(t193) * t358 - t147 * sin(t193) + t202;
t6 = t18 * t332 - t75 * t12 - t116 * t15 - 0.4e1 * ((t66 * t302 / 0.2e1 + t60 * t356) * t183 + t329 * t357 + (t161 - 0.1e1 / 0.2e1) * t66 * t147) * t60 + (-t119 * t381 - t239 * t144) * t184 - (t239 * t119 - t144 * t381) * t178;
t199 = t381 * t178 - t239 * t184;
t9 = t18 * t365 - t116 * t12 - t141 * t15 + t63 * (t161 * t369 + t183 * t302 - t147) + (((t240 * t164 - t308) * rSges(3,1) + t199 * rSges(3,2)) * t183 + t177 * (t300 + (t125 * t293 + t308) * rSges(3,2) + t199 * rSges(3,1))) * m(3);
t230 = -t9 * t326 - t93 * t6;
t220 = t82 * t316 + t73 * t91;
t129 = pkin(5) * t174 + t180 * t350;
t217 = pkin(2) * t292 - t126 * t166;
t76 = t165 * t129 + t217 * t163;
t67 = t108 * t150 + t153 * t76;
t28 = -t220 * t255 + t67 * t259;
t214 = t114 * t91 + t82 * t301;
t40 = -t214 * t255 + t67 * t271;
t229 = t28 * t91 + t40 * t328;
t68 = t108 * t153 - t150 * t76;
t29 = t220 * t256 + t68 * t259;
t41 = t214 * t256 + t68 * t271;
t228 = t29 * t91 + t41 * t328;
t219 = t83 * t315 + t74 * t92;
t130 = pkin(5) * t176 + t182 * t349;
t216 = pkin(2) * t291 - t127 * t166;
t77 = t165 * t130 + t216 * t163;
t69 = t109 * t151 + t154 * t77;
t30 = -t219 * t253 + t69 * t258;
t213 = t115 * t92 + t83 * t301;
t42 = -t213 * t253 + t69 * t270;
t227 = t30 * t92 + t42 * t327;
t70 = t109 * t154 - t151 * t77;
t31 = t219 * t254 + t70 * t258;
t43 = t213 * t254 + t70 * t270;
t226 = t31 * t92 + t43 * t327;
t218 = t84 * t314 + t75 * t93;
t131 = pkin(5) * t178 + t184 * t348;
t215 = pkin(2) * t290 - t128 * t166;
t78 = t165 * t131 + t215 * t163;
t71 = t110 * t152 + t155 * t78;
t32 = -t218 * t251 + t71 * t257;
t212 = t116 * t93 + t84 * t301;
t44 = -t212 * t251 + t71 * t269;
t225 = t32 * t93 + t44 * t326;
t72 = t110 * t155 - t152 * t78;
t33 = t218 * t252 + t72 * t257;
t45 = t212 * t252 + t72 * t269;
t224 = t45 * t326 + t33 * t93;
t208 = t82 * t268 + t91 * t334;
t207 = t83 * t267 + t92 * t333;
t206 = t84 * t266 + t93 * t332;
t170 = xDDP(3);
t79 = -t163 * t129 + t217 * t165;
t205 = t170 * t79 + t171 * t68 + t172 * t67;
t80 = -t163 * t130 + t216 * t165;
t204 = t170 * t80 + t171 * t70 + t172 * t69;
t81 = -t163 * t131 + t215 * t165;
t203 = t170 * t81 + t171 * t72 + t172 * t71;
t57 = t60 ^ 2;
t56 = t59 ^ 2;
t55 = t58 ^ 2;
t54 = t81 * t269 + (t116 * t90 + t87 * t301) * t320;
t53 = t80 * t270 + (t115 * t89 + t86 * t301) * t321;
t52 = t79 * t271 + (t114 * t88 + t85 * t301) * t322;
t48 = t81 * t257 + (t87 * t314 + t75 * t90) * t320;
t47 = t80 * t258 + (t86 * t315 + t74 * t89) * t321;
t46 = t79 * t259 + (t85 * t316 + t73 * t88) * t322;
t3 = (-t12 * t99 + (-t144 * t184 - t178 * t187) * t63) * t164 + (t18 - t122) * t156 + (-t96 * t15 + (-0.2e1 * t184 * (rSges(3,1) * t329 + t60 * t336) * t66 + t375 * t178 * (t63 + t57)) * t164 - t57 * t166 * t134) * m(3);
t2 = (-t11 * t98 + (-t144 * t182 - t176 * t187) * t62) * t164 + (-t17 - t121) * t156 + (-t95 * t14 + (-0.2e1 * t182 * (rSges(3,1) * t330 + t59 * t337) * t65 + t374 * t176 * (t62 + t56)) * t164 - t56 * t166 * t133) * m(3);
t1 = (-t10 * t97 + (-t144 * t180 - t174 * t187) * t61) * t164 + (-t16 - t120) * t156 + (-t94 * t13 + (-0.2e1 * t180 * (rSges(3,1) * t331 + t58 * t338) * t64 + t373 * t174 * (t61 + t55)) * t164 - t55 * t166 * t132) * m(3);
t19 = [(-g(1) + t172) * m(4) + (t71 * t3 + t203 * (-t206 * t251 + t71 * t317)) * t107 + (t69 * t2 + t204 * (-t207 * t253 + t69 * t318)) * t370 + (t67 * t1 + t205 * (-t208 * t255 + t67 * t319)) * t371 + ((t32 * t90 + t44 * t323) * t170 + t225 * t297 + (-t225 * t172 + t230) * t155) * t320 + ((t30 * t89 + t42 * t324) * t170 + t227 * t298 + (-t227 * t172 + t231) * t154) * t321 + ((t28 * t88 + t40 * t325) * t170 + t229 * t299 + (-t229 * t172 + t232) * t153) * t322; (-g(2) + t171) * m(4) + (t72 * t3 + t203 * (t206 * t252 + t72 * t317)) * t107 + (t70 * t2 + t204 * (t207 * t254 + t70 * t318)) * t370 + (t68 * t1 + t205 * (t208 * t256 + t68 * t319)) * t371 + ((t45 * t323 + t33 * t90) * t170 - t224 * t294 + (t224 * t171 - t230) * t152) * t320 + ((t31 * t89 + t43 * t324) * t170 - t226 * t295 + (t226 * t171 - t231) * t151) * t321 + ((t29 * t88 + t41 * t325) * t170 - t228 * t296 + (t228 * t171 - t232) * t150) * t322; (-g(3) + t170) * m(4) + (t81 * t3 + t203 * (t81 * t317 + (t87 * t266 + t90 * t332) * t320)) * t107 + (t80 * t2 + t204 * (t80 * t318 + (t86 * t267 + t89 * t333) * t321)) * t370 + (t79 * t1 + t205 * (t79 * t319 + (t85 * t268 + t88 * t334) * t322)) * t371 + ((t54 * t323 + t48 * t90) * t170 + t90 * t6 + t9 * t323 + (-t294 + t297) * (t54 * t326 + t48 * t93)) * t320 + ((t53 * t324 + t47 * t89) * t170 + t89 * t5 + t8 * t324 + (-t295 + t298) * (t53 * t327 + t47 * t92)) * t321 + ((t52 * t325 + t46 * t88) * t170 + t88 * t4 + t7 * t325 + (-t296 + t299) * (t52 * t328 + t46 * t91)) * t322;];
tauX  = t19;
