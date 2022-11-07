% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR8V2G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2022-11-07 13:12
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:10:50
% EndTime: 2022-11-07 13:10:57
% DurationCPUTime: 7.27s
% Computational Cost: add. (30216->580), mult. (34284->906), div. (5283->9), fcn. (22491->62), ass. (0->385)
t447 = 2 * pkin(1);
t390 = pkin(5) + qJ(3,3);
t182 = -pkin(6) - t390;
t171 = 0.1e1 / t182;
t215 = cos(qJ(2,3));
t175 = t215 * pkin(2);
t189 = qJ(2,3) + pkin(7);
t158 = cos(t189);
t401 = pkin(3) * t158;
t107 = t175 + t401;
t100 = 0.1e1 / t107;
t227 = xDP(3);
t373 = t100 * t227;
t148 = sin(t189);
t209 = sin(qJ(2,3));
t405 = pkin(2) * t209;
t110 = pkin(3) * t148 + t405;
t234 = 0.2e1 * qJ(2,3);
t188 = t234 + pkin(7);
t147 = sin(t188);
t432 = pkin(2) * pkin(3);
t336 = 0.2e1 * t432;
t194 = sin(t234);
t246 = pkin(2) ^ 2;
t339 = t246 * t194;
t134 = 0.2e1 * t189;
t123 = sin(t134);
t244 = pkin(3) ^ 2;
t341 = t244 * t123;
t263 = t147 * t336 + t339 + t341;
t70 = t110 * t447 + t263;
t313 = t70 * t373;
t229 = xDP(1);
t203 = legFrame(3,3);
t165 = sin(t203);
t168 = cos(t203);
t139 = t175 + pkin(1);
t210 = sin(qJ(1,3));
t216 = cos(qJ(1,3));
t78 = t139 * t210 + t182 * t216;
t81 = t139 * t216 - t182 * t210;
t85 = -t165 * t210 + t168 * t216;
t59 = -t165 * t78 + t168 * t81 + t401 * t85;
t380 = t59 * t229;
t228 = xDP(2);
t84 = t165 * t216 + t168 * t210;
t58 = t165 * t81 + t168 * t78 + t401 * t84;
t381 = t58 * t228;
t31 = (t380 + t381 + t313 / 0.2e1) * t171;
t441 = -0.2e1 * t31;
t391 = pkin(5) + qJ(3,2);
t183 = -pkin(6) - t391;
t172 = 0.1e1 / t183;
t217 = cos(qJ(2,2));
t176 = t217 * pkin(2);
t191 = qJ(2,2) + pkin(7);
t161 = cos(t191);
t400 = pkin(3) * t161;
t108 = t176 + t400;
t102 = 0.1e1 / t108;
t370 = t102 * t227;
t150 = sin(t191);
t211 = sin(qJ(2,2));
t404 = pkin(2) * t211;
t111 = pkin(3) * t150 + t404;
t236 = 0.2e1 * qJ(2,2);
t190 = t236 + pkin(7);
t149 = sin(t190);
t195 = sin(t236);
t338 = t246 * t195;
t135 = 0.2e1 * t191;
t124 = sin(t135);
t340 = t244 * t124;
t262 = t149 * t336 + t338 + t340;
t71 = t111 * t447 + t262;
t312 = t71 * t370;
t204 = legFrame(2,3);
t166 = sin(t204);
t169 = cos(t204);
t140 = t176 + pkin(1);
t212 = sin(qJ(1,2));
t218 = cos(qJ(1,2));
t79 = t140 * t212 + t183 * t218;
t82 = t140 * t218 - t183 * t212;
t87 = -t166 * t212 + t169 * t218;
t61 = -t166 * t79 + t169 * t82 + t400 * t87;
t378 = t61 * t229;
t86 = t166 * t218 + t169 * t212;
t60 = t166 * t82 + t169 * t79 + t400 * t86;
t379 = t60 * t228;
t32 = (t378 + t379 + t312 / 0.2e1) * t172;
t440 = -0.2e1 * t32;
t389 = qJ(3,1) + pkin(5);
t184 = -pkin(6) - t389;
t173 = 0.1e1 / t184;
t219 = cos(qJ(2,1));
t177 = t219 * pkin(2);
t193 = qJ(2,1) + pkin(7);
t163 = cos(t193);
t399 = pkin(3) * t163;
t109 = t177 + t399;
t104 = 0.1e1 / t109;
t367 = t104 * t227;
t152 = sin(t193);
t213 = sin(qJ(2,1));
t403 = pkin(2) * t213;
t112 = pkin(3) * t152 + t403;
t238 = 0.2e1 * qJ(2,1);
t192 = t238 + pkin(7);
t151 = sin(t192);
t196 = sin(t238);
t342 = t196 * t246;
t136 = 0.2e1 * t193;
t125 = sin(t136);
t356 = t125 * t244;
t261 = t151 * t336 + t342 + t356;
t72 = t112 * t447 + t261;
t311 = t72 * t367;
t205 = legFrame(1,3);
t167 = sin(t205);
t170 = cos(t205);
t141 = t177 + pkin(1);
t214 = sin(qJ(1,1));
t220 = cos(qJ(1,1));
t80 = t141 * t214 + t184 * t220;
t83 = t141 * t220 - t184 * t214;
t89 = -t167 * t214 + t170 * t220;
t63 = -t167 * t80 + t170 * t83 + t399 * t89;
t376 = t63 * t229;
t88 = t167 * t220 + t170 * t214;
t62 = t167 * t83 + t170 * t80 + t399 * t88;
t377 = t62 * t228;
t33 = (t376 + t377 + t311 / 0.2e1) * t173;
t439 = -0.2e1 * t33;
t162 = cos(t192);
t446 = rSges(3,1) * t151 + rSges(3,2) * t162;
t160 = cos(t190);
t445 = rSges(3,1) * t149 + rSges(3,2) * t160;
t157 = cos(t188);
t444 = rSges(3,1) * t147 + rSges(3,2) * t157;
t443 = -2 * pkin(1);
t202 = cos(pkin(7));
t318 = t202 * t432;
t413 = t246 / 0.2e1;
t438 = -0.4e1 * pkin(1) * (t318 + t413 + t244 / 0.2e1);
t200 = t227 ^ 2;
t437 = 0.1e1 / t107 ^ 2;
t436 = 0.1e1 / t108 ^ 2;
t435 = 0.1e1 / t109 ^ 2;
t245 = pkin(2) * t246;
t402 = pkin(2) * t244;
t434 = -0.2e1 * t245 - 0.4e1 * t402;
t433 = m(3) * pkin(1);
t230 = pkin(2) * m(3);
t363 = t110 * t100;
t55 = (t227 * t363 + t228 * t84 + t229 * t85) * t171;
t431 = -t55 / 0.4e1;
t362 = t111 * t102;
t56 = (t227 * t362 + t228 * t86 + t229 * t87) * t172;
t430 = -t56 / 0.4e1;
t361 = t112 * t104;
t57 = (t227 * t361 + t228 * t88 + t229 * t89) * t173;
t429 = -t57 / 0.4e1;
t428 = t70 / 0.2e1;
t427 = t71 / 0.2e1;
t426 = t72 / 0.2e1;
t226 = m(2) * rSges(2,1);
t425 = m(2) * rSges(2,2);
t424 = m(3) * rSges(3,2);
t74 = -rSges(3,1) * t158 + rSges(3,2) * t148 - t139;
t423 = m(3) * t74;
t75 = -rSges(3,1) * t161 + rSges(3,2) * t150 - t140;
t422 = m(3) * t75;
t76 = -rSges(3,1) * t163 + rSges(3,2) * t152 - t141;
t421 = m(3) * t76;
t420 = rSges(2,2) * pkin(1);
t137 = t226 + t230;
t201 = sin(pkin(7));
t317 = t201 * t424;
t385 = rSges(3,1) * t202;
t77 = -m(3) * t385 - t137 + t317;
t419 = g(3) * t77;
t240 = rSges(2,2) ^ 2;
t242 = (rSges(2,1) ^ 2);
t106 = t246 * m(3) + ((-t240 + t242) * m(2)) + Icges(2,2) - Icges(2,1);
t418 = t106 / 0.2e1;
t239 = rSges(3,2) ^ 2;
t241 = rSges(3,1) ^ 2;
t115 = m(3) * (-t239 + t241) - Icges(3,1) + Icges(3,2);
t417 = t115 / 0.2e1;
t416 = -t227 / 0.2e1;
t414 = -0.3e1 / 0.4e1 * t246;
t221 = (rSges(2,3) + pkin(5));
t412 = m(2) * t221;
t411 = m(3) * t171;
t410 = m(3) * t172;
t409 = m(3) * t173;
t179 = (rSges(3,3) + t390);
t408 = m(3) * t179;
t180 = (rSges(3,3) + t391);
t407 = m(3) * t180;
t181 = (rSges(3,3) + t389);
t406 = m(3) * t181;
t398 = pkin(3) * t246;
t101 = t100 * t437;
t126 = cos(t134);
t178 = t244 + t246;
t197 = cos(t234);
t270 = t126 * t244 + t197 * t246 + t178;
t301 = -0.2e1 * t318;
t337 = -0.2e1 * t432;
t360 = (0.2e1 * t318 + t178) * t200;
t375 = t100 * t171;
t52 = pkin(1) * t55;
t10 = t101 * t171 * t360 + (-(t157 * t337 - t270 + t301) * t55 / 0.2e1 + (t441 + t52) * t107) * t55 * t375;
t117 = rSges(3,2) * t408 - Icges(3,6);
t120 = rSges(3,1) * t408 - Icges(3,5);
t300 = m(1) * rSges(1,1) + m(2) * pkin(1) + t433;
t92 = t425 + (rSges(3,1) * t201 + rSges(3,2) * t202) * m(3);
t266 = -t209 * t92 - t215 * t77 + t300;
t316 = t55 * t373;
t291 = -0.2e1 * t316;
t129 = m(1) * rSges(1,2) - t412;
t296 = t129 - t408;
t306 = t101 * t200 * t110;
t308 = t424 * t447;
t309 = -0.2e1 * rSges(3,1) * t433;
t310 = 2 * m(2) * t420;
t321 = pkin(2) * t408;
t327 = t137 * t443;
t335 = rSges(2,1) * t412 - Icges(2,5);
t142 = -rSges(3,1) * t424 + Icges(3,4);
t143 = rSges(2,1) * t425 - Icges(2,4);
t133 = pkin(2) * t317;
t256 = Icges(1,3) + ((rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1)) + Icges(2,1) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(2,2) / 0.2e1 + Icges(3,2) / 0.2e1 - t133;
t247 = pkin(1) ^ 2;
t257 = (pkin(5) ^ 2) + ((2 * pkin(5) + rSges(2,3)) * rSges(2,3)) + t247 + t240 / 0.2e1 + t242 / 0.2e1;
t267 = t247 + t239 / 0.2e1 + t241 / 0.2e1 + t413;
t326 = t215 * t447;
t329 = -2 * t420;
t334 = t157 + t202;
t34 = t126 * t417 + t197 * t418 + t137 * t326 + (t209 * t329 + t257) * m(2) + t142 * t123 - t143 * t194 + ((t179 ^ 2) + (-pkin(2) * t147 + t148 * t443) * rSges(3,2) + (pkin(2) * t334 + t158 * t447) * rSges(3,1) + t267) * m(3) + t256;
t351 = t143 * t197;
t354 = t142 * t126;
t355 = (rSges(2,2) * t412 - Icges(2,6)) * t227;
t359 = t115 * t123;
t366 = t106 * t194;
t271 = t221 * t425 - Icges(2,6);
t272 = -t221 * t226 + Icges(2,5);
t49 = (t117 * t201 - t120 * t202 + t272 - t321) * t209 - (t117 * t202 + t120 * t201 + t271) * t215;
t232 = 0.2e1 * pkin(7);
t153 = cos(t232 + qJ(2,3));
t159 = cos(qJ(2,3) - pkin(7));
t22 = -t52 - (-t380 / 0.2e1 - t381 / 0.2e1 - t313 / 0.4e1) * t171;
t233 = 0.3e1 * qJ(2,3);
t243 = pkin(3) * t244;
t285 = 0.3e1 / 0.4e1 * t246;
t286 = 0.3e1 / 0.4e1 * t244;
t293 = -0.2e1 * t243 - 0.4e1 * t398;
t299 = -t375 / 0.2e1;
t322 = -0.2e1 * t398;
t323 = -0.2e1 * t402;
t331 = -0.6e1 * t246 - 0.3e1 * t244;
t46 = (t182 ^ 2 + t247) * t55;
t64 = t182 * t373;
t7 = (t153 * t323 + t158 * t293 + t159 * t322 + t215 * t434 + t438) * t200 * t437 * t299 + (t263 * t437 * t182 * t416 + ((-t243 * cos(0.3e1 * t189) - t245 * cos(t233)) * t431 - (t341 / 0.2e1 + t339 / 0.2e1) * t64 + (-(-cos(t233 + pkin(7)) - t159) * t55 * t285 + (pkin(1) * t441 + t331 * t431 + t46) * t158) * pkin(3) + ((-t414 * t55 + t46) * t215 - (-cos(t232 + t233) - t153 - 0.2e1 * t215) * t55 * t286 + (-t64 * t147 - 0.2e1 * t22 * t334) * pkin(3) - (pkin(3) * t334 + t326) * t31) * pkin(2) + (-t31 / 0.2e1 - t22) * t270) * t100) * t171 * t55;
t93 = -g(1) * t165 + g(2) * t168;
t96 = g(1) * t168 + g(2) * t165;
t397 = (-t34 * t10 + t49 * t306 - t7 * t423 + t291 * t354 - t55 * t408 * t441 + (-t266 * t93 + t296 * t96) * t216 + (t266 * t96 + t296 * t93) * t210 - (-t359 - t366) * t316 + (t100 * t355 * t209 + (-t120 * t158 + t117 * t148 - (t321 + t335) * t215) * t373 - (t309 * t148 - t308 * t158 + t327 * t209 - t310 * t215) * t55) * t373 + 0.2e1 * (t444 * t230 + t351) * t316) * t171;
t278 = rSges(3,1) * t148 + rSges(3,2) * t158;
t4 = (-t55 ^ 2 * t179 - t74 * t10 - t210 * t96 + t216 * t93 - t7 + (t278 + t405) * t291) * m(3);
t396 = t171 * t4;
t103 = t102 * t436;
t127 = cos(t135);
t198 = cos(t236);
t269 = t127 * t244 + t198 * t246 + t178;
t372 = t102 * t172;
t53 = pkin(1) * t56;
t11 = t103 * t172 * t360 + (-(t160 * t337 - t269 + t301) * t56 / 0.2e1 + (t440 + t53) * t108) * t56 * t372;
t118 = rSges(3,2) * t407 - Icges(3,6);
t121 = rSges(3,1) * t407 - Icges(3,5);
t265 = -t211 * t92 - t217 * t77 + t300;
t315 = t56 * t370;
t289 = -0.2e1 * t315;
t295 = t129 - t407;
t304 = t103 * t200 * t111;
t320 = pkin(2) * t407;
t325 = t217 * t447;
t333 = t160 + t202;
t35 = t127 * t417 + t198 * t418 + t137 * t325 + (t211 * t329 + t257) * m(2) + t142 * t124 - t143 * t195 + ((t180 ^ 2) + (-pkin(2) * t149 + t150 * t443) * rSges(3,2) + (pkin(2) * t333 + t161 * t447) * rSges(3,1) + t267) * m(3) + t256;
t350 = t143 * t198;
t353 = t142 * t127;
t358 = t115 * t124;
t365 = t106 * t195;
t50 = (t118 * t201 - t121 * t202 + t272 - t320) * t211 - (t118 * t202 + t121 * t201 + t271) * t217;
t154 = cos(t232 + qJ(2,2));
t156 = cos(-pkin(7) + qJ(2,2));
t23 = -t53 - (-t378 / 0.2e1 - t379 / 0.2e1 - t312 / 0.4e1) * t172;
t235 = 0.3e1 * qJ(2,2);
t298 = -t372 / 0.2e1;
t47 = (t183 ^ 2 + t247) * t56;
t65 = t183 * t370;
t8 = (t154 * t323 + t156 * t322 + t161 * t293 + t217 * t434 + t438) * t200 * t436 * t298 + (t262 * t436 * t183 * t416 + ((-t243 * cos(0.3e1 * t191) - t245 * cos(t235)) * t430 - (t340 / 0.2e1 + t338 / 0.2e1) * t65 + (-(-cos(pkin(7) + t235) - t156) * t56 * t285 + (pkin(1) * t440 + t331 * t430 + t47) * t161) * pkin(3) + ((-t414 * t56 + t47) * t217 - (-cos(t232 + t235) - t154 - 0.2e1 * t217) * t56 * t286 + (-t65 * t149 - 0.2e1 * t23 * t333) * pkin(3) - (pkin(3) * t333 + t325) * t32) * pkin(2) + (-t32 / 0.2e1 - t23) * t269) * t102) * t172 * t56;
t94 = -g(1) * t166 + g(2) * t169;
t97 = g(1) * t169 + g(2) * t166;
t395 = t172 * (-t35 * t11 + t50 * t304 - t8 * t422 + t289 * t353 - t56 * t407 * t440 + (-t265 * t94 + t295 * t97) * t218 + (t265 * t97 + t295 * t94) * t212 - (-t358 - t365) * t315 + (t102 * t355 * t211 + (-t121 * t161 + t118 * t150 - (t320 + t335) * t217) * t370 - (t309 * t150 - t308 * t161 + t327 * t211 - t310 * t217) * t56) * t370 + 0.2e1 * (t445 * t230 + t350) * t315);
t277 = rSges(3,1) * t150 + rSges(3,2) * t161;
t5 = (-t56 ^ 2 * t180 - t75 * t11 - t212 * t97 + t218 * t94 - t8 + (t277 + t404) * t289) * m(3);
t394 = t172 * t5;
t119 = rSges(3,2) * t406 - Icges(3,6);
t105 = t104 * t435;
t128 = cos(t136);
t199 = cos(t238);
t268 = t128 * t244 + t199 * t246 + t178;
t369 = t104 * t173;
t54 = pkin(1) * t57;
t12 = t105 * t173 * t360 + (-(t162 * t337 - t268 + t301) * t57 / 0.2e1 + (t439 + t54) * t109) * t57 * t369;
t122 = rSges(3,1) * t406 - Icges(3,5);
t264 = -t213 * t92 - t219 * t77 + t300;
t314 = t57 * t367;
t287 = -0.2e1 * t314;
t294 = t129 - t406;
t302 = t105 * t200 * t112;
t319 = pkin(2) * t406;
t349 = t143 * t199;
t352 = t142 * t128;
t357 = t115 * t125;
t324 = t219 * t447;
t332 = t162 + t202;
t36 = t128 * t417 + t199 * t418 + t137 * t324 + (t213 * t329 + t257) * m(2) + t142 * t125 - t143 * t196 + ((t181 ^ 2) + (-pkin(2) * t151 + t152 * t443) * rSges(3,2) + (pkin(2) * t332 + t163 * t447) * rSges(3,1) + t267) * m(3) + t256;
t364 = t106 * t196;
t51 = (t119 * t201 - t122 * t202 + t272 - t319) * t213 - (t119 * t202 + t122 * t201 + t271) * t219;
t155 = cos(t232 + qJ(2,1));
t164 = cos(qJ(2,1) - pkin(7));
t237 = 0.3e1 * qJ(2,1);
t24 = -t54 - (-t376 / 0.2e1 - t377 / 0.2e1 - t311 / 0.4e1) * t173;
t297 = -t369 / 0.2e1;
t48 = (t184 ^ 2 + t247) * t57;
t66 = t184 * t367;
t9 = (t155 * t323 + t163 * t293 + t164 * t322 + t219 * t434 + t438) * t200 * t435 * t297 + (t261 * t435 * t184 * t416 + ((-t243 * cos(0.3e1 * t193) - t245 * cos(t237)) * t429 - (t356 / 0.2e1 + t342 / 0.2e1) * t66 + (-(-cos(t237 + pkin(7)) - t164) * t57 * t285 + (pkin(1) * t439 + t331 * t429 + t48) * t163) * pkin(3) + ((-t414 * t57 + t48) * t219 - (-cos(t237 + t232) - t155 - 0.2e1 * t219) * t57 * t286 + (-t66 * t151 - 0.2e1 * t24 * t332) * pkin(3) - (pkin(3) * t332 + t324) * t33) * pkin(2) + (-t33 / 0.2e1 - t24) * t268) * t104) * t173 * t57;
t95 = -g(1) * t167 + g(2) * t170;
t98 = g(1) * t170 + g(2) * t167;
t393 = t173 * (-t36 * t12 + t51 * t302 - t9 * t421 + t287 * t352 - t57 * t406 * t439 + (-t264 * t95 + t294 * t98) * t220 + (t264 * t98 + t294 * t95) * t214 - (-t357 - t364) * t314 + (t104 * t355 * t213 + (-t122 * t163 + t119 * t152 - (t319 + t335) * t219) * t367 - (t309 * t152 - t308 * t163 + t327 * t213 - t310 * t219) * t57) * t367 + 0.2e1 * (t446 * t230 + t349) * t314);
t276 = rSges(3,1) * t152 + rSges(3,2) * t163;
t6 = (-t57 ^ 2 * t181 - t76 * t12 - t214 * t98 + t220 * t95 - t9 + (t276 + t403) * t287) * m(3);
t392 = t173 * t6;
t206 = xDDP(3);
t374 = t100 * t206;
t371 = t102 * t206;
t368 = t104 * t206;
t207 = xDDP(2);
t348 = t171 * t207;
t208 = xDDP(1);
t347 = t171 * t208;
t346 = t172 * t207;
t345 = t172 * t208;
t344 = t173 * t207;
t343 = t173 * t208;
t307 = t171 * t374;
t305 = t172 * t371;
t303 = t173 * t368;
t275 = t210 * t93 + t216 * t96;
t274 = t212 * t94 + t218 * t97;
t273 = t214 * t95 + t220 * t98;
t90 = g(3) * t92;
t73 = -0.2e1 * t133 + ((t240 + t242) * m(2)) + Icges(2,3) + Icges(3,3) + (0.2e1 * pkin(2) * t385 + t239 + t241 + t246) * m(3);
t45 = (t112 * t76 + t426) * m(3) * t369;
t44 = (t111 * t75 + t427) * m(3) * t372;
t43 = (t110 * t74 + t428) * m(3) * t375;
t42 = (t76 * t89 + t63) * t409;
t41 = (t76 * t88 + t62) * t409;
t40 = (t75 * t87 + t61) * t410;
t39 = (t75 * t86 + t60) * t410;
t38 = (t74 * t85 + t59) * t411;
t37 = (t74 * t84 + t58) * t411;
t21 = (t36 * t89 + t421 * t63) * t173;
t20 = (t36 * t88 + t421 * t62) * t173;
t19 = (t35 * t87 + t422 * t61) * t172;
t18 = (t35 * t86 + t422 * t60) * t172;
t17 = (t34 * t85 + t423 * t59) * t171;
t16 = (t34 * t84 + t423 * t58) * t171;
t15 = (t51 - (t112 * t36 + t421 * t426) * t173) * t104;
t14 = (t50 - (t111 * t35 + t422 * t427) * t172) * t102;
t13 = (t49 - (t110 * t34 + t423 * t428) * t171) * t100;
t1 = [-(-t21 * t89 - t42 * t63) * t343 - (-t21 * t88 - t42 * t62) * t344 - (-t21 * t112 - t42 * t426 + t89 * t51) * t303 - t89 * t393 - t63 * t392 - (-t19 * t87 - t40 * t61) * t345 - (-t19 * t86 - t40 * t60) * t346 - (-t19 * t111 - t40 * t427 + t87 * t50) * t305 - t87 * t395 - t61 * t394 - (-t17 * t85 - t38 * t59) * t347 - (-t17 * t84 - t38 * t58) * t348 - (-t17 * t110 - t38 * t428 + t85 * t49) * t307 - t85 * t397 - t59 * t396 + (-g(1) + t208) * m(4); -(-t20 * t89 - t41 * t63) * t343 - (-t20 * t88 - t41 * t62) * t344 - (-t20 * t112 - t41 * t426 + t88 * t51) * t303 - t88 * t393 - t62 * t392 - (-t18 * t87 - t39 * t61) * t345 - (-t18 * t86 - t39 * t60) * t346 - (-t18 * t111 - t39 * t427 + t86 * t50) * t305 - t86 * t395 - t60 * t394 - (-t16 * t85 - t37 * t59) * t347 - (-t16 * t84 - t37 * t58) * t348 - (-t16 * t110 - t37 * t428 + t84 * t49) * t307 - t84 * t397 - t58 * t396 + (-g(2) + t207) * m(4); -(t15 * t89 - t45 * t63) * t343 - (t15 * t88 - t45 * t62) * t344 + (t104 * t73 - (-t45 * t426 + (t104 * t51 + t15) * t112) * t173) * t368 - t361 * t393 + t104 * (-t51 * t12 + t73 * t302 + (-t273 * t77 + t90) * t213 + t219 * (t273 * t92 + t419) - (-(t357 / 0.2e1 - t352 + t364 / 0.2e1 + t349 + (t137 * t213 + t219 * t425) * pkin(1)) * t57 + (t276 * (-t54 - (-t311 - 0.2e1 * t376 - 0.2e1 * t377) * t173) + (0.2e1 * t33 * t213 - t446 * t57) * pkin(2)) * m(3)) * t57) + t72 * t6 * t297 - (t14 * t87 - t44 * t61) * t345 - (t14 * t86 - t44 * t60) * t346 + (t102 * t73 - (-t44 * t427 + (t102 * t50 + t14) * t111) * t172) * t371 - t362 * t395 + t102 * (-t50 * t11 + t73 * t304 + (-t274 * t77 + t90) * t211 + t217 * (t274 * t92 + t419) - (-(t358 / 0.2e1 - t353 + t365 / 0.2e1 + t350 + (t137 * t211 + t217 * t425) * pkin(1)) * t56 + (t277 * (-t53 - (-t312 - 0.2e1 * t378 - 0.2e1 * t379) * t172) + (0.2e1 * t32 * t211 - t445 * t56) * pkin(2)) * m(3)) * t56) + t71 * t5 * t298 - (t13 * t85 - t43 * t59) * t347 - (t13 * t84 - t43 * t58) * t348 + (t100 * t73 - (-t43 * t428 + (t100 * t49 + t13) * t110) * t171) * t374 - t363 * t397 + t100 * (-t49 * t10 + t73 * t306 + (-t275 * t77 + t90) * t209 + t215 * (t275 * t92 + t419) - (-(t359 / 0.2e1 - t354 + t366 / 0.2e1 + t351 + (t137 * t209 + t215 * t425) * pkin(1)) * t55 + (t278 * (-t52 - (-t313 - 0.2e1 * t380 - 0.2e1 * t381) * t171) + (0.2e1 * t31 * t209 - t444 * t55) * pkin(2)) * m(3)) * t55) + t70 * t4 * t299 + (-g(3) + t206) * m(4);];
tauX  = t1;
