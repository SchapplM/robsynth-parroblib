% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRR1G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
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
%   pkin=[a2,a3,d2,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRR1G1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1G1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:35
% EndTime: 2019-05-03 14:47:35
% DurationCPUTime: 0.60s
% Computational Cost: add. (694->124), mult. (1360->271), div. (462->4), fcn. (1318->14), ass. (0->124)
t294 = xDP(3);
t284 = t294 ^ 2;
t288 = sin(qJ(2,3));
t291 = cos(qJ(2,3));
t269 = t291 * rSges(2,1) - t288 * rSges(2,2);
t343 = m(2) * t269;
t289 = sin(qJ(2,2));
t292 = cos(qJ(2,2));
t270 = t292 * rSges(2,1) - t289 * rSges(2,2);
t342 = m(2) * t270;
t290 = sin(qJ(2,1));
t293 = cos(qJ(2,1));
t271 = t293 * rSges(2,1) - t290 * rSges(2,2);
t341 = m(2) * t271;
t307 = 0.1e1 / pkin(2);
t340 = m(2) * t307;
t297 = m(1) + m(2);
t339 = pkin(2) * t297;
t338 = m(3) * t284;
t298 = xP(3);
t279 = sin(t298);
t280 = cos(t298);
t301 = koppelP(3,2);
t304 = koppelP(3,1);
t260 = t279 * t304 + t280 * t301;
t263 = -t279 * t301 + t280 * t304;
t285 = legFrame(3,3);
t273 = sin(t285);
t276 = cos(t285);
t295 = xDP(2);
t296 = xDP(1);
t281 = 0.1e1 / t288;
t326 = t281 * t307;
t236 = (-t276 * (-t260 * t294 + t296) - t273 * (t263 * t294 + t295)) * t326;
t233 = t236 ^ 2;
t337 = t233 * t281;
t302 = koppelP(2,2);
t305 = koppelP(2,1);
t261 = t279 * t305 + t280 * t302;
t264 = -t279 * t302 + t280 * t305;
t286 = legFrame(2,3);
t274 = sin(t286);
t277 = cos(t286);
t282 = 0.1e1 / t289;
t323 = t282 * t307;
t237 = (-t277 * (-t261 * t294 + t296) - t274 * (t264 * t294 + t295)) * t323;
t234 = t237 ^ 2;
t336 = t234 * t282;
t303 = koppelP(1,2);
t306 = koppelP(1,1);
t262 = t279 * t306 + t280 * t303;
t265 = -t279 * t303 + t280 * t306;
t287 = legFrame(1,3);
t275 = sin(t287);
t278 = cos(t287);
t283 = 0.1e1 / t290;
t320 = t283 * t307;
t238 = (-t278 * (-t262 * t294 + t296) - t275 * (t265 * t294 + t295)) * t320;
t235 = t238 ^ 2;
t335 = t235 * t283;
t334 = t273 * t307;
t333 = t274 * t307;
t332 = t275 * t307;
t331 = t276 * t307;
t330 = t277 * t307;
t329 = t278 * t307;
t221 = (-t291 * t343 + t339) * t337;
t328 = t281 * t221;
t327 = t281 * t284;
t222 = (-t292 * t342 + t339) * t336;
t325 = t282 * t222;
t324 = t282 * t284;
t223 = (-t293 * t341 + t339) * t335;
t322 = t283 * t223;
t321 = t283 * t284;
t319 = m(2) * t233 * (t288 * rSges(2,1) + t291 * rSges(2,2));
t318 = m(2) * t234 * (t289 * rSges(2,1) + t292 * rSges(2,2));
t317 = m(2) * t235 * (t290 * rSges(2,1) + t293 * rSges(2,2));
t316 = t269 * t340;
t315 = t270 * t340;
t314 = t271 * t340;
t272 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(2,3);
t218 = (pkin(2) * t343 - t272 * t291) * t337;
t313 = t218 * t326;
t219 = (pkin(2) * t342 - t272 * t292) * t336;
t312 = t219 * t323;
t220 = (pkin(2) * t341 - t272 * t293) * t335;
t311 = t220 * t320;
t310 = t281 * t319;
t309 = t282 * t318;
t308 = t283 * t317;
t300 = rSges(3,1);
t299 = rSges(3,2);
t259 = t275 * t293 + t290 * t278;
t258 = -t275 * t290 + t278 * t293;
t257 = t274 * t292 + t289 * t277;
t256 = -t274 * t289 + t277 * t292;
t255 = t273 * t291 + t288 * t276;
t254 = -t273 * t288 + t276 * t291;
t253 = (t259 * t297 - t275 * t314) * t283;
t252 = (t258 * t297 - t278 * t314) * t283;
t251 = (t257 * t297 - t274 * t315) * t282;
t250 = (t256 * t297 - t277 * t315) * t282;
t249 = (t255 * t297 - t273 * t316) * t281;
t248 = (t254 * t297 - t276 * t316) * t281;
t247 = (t278 * t262 - t275 * t265) * t320;
t246 = (t277 * t261 - t274 * t264) * t323;
t245 = (t276 * t260 - t273 * t263) * t326;
t244 = (t259 * t341 - t272 * t332) * t283;
t243 = (t258 * t341 - t272 * t329) * t283;
t242 = (t257 * t342 - t272 * t333) * t282;
t241 = (t256 * t342 - t272 * t330) * t282;
t240 = (t255 * t343 - t272 * t334) * t281;
t239 = (t254 * t343 - t272 * t331) * t281;
t232 = (-t258 * t262 + t259 * t265) * t283;
t231 = (-t256 * t261 + t257 * t264) * t282;
t230 = (-t254 * t260 + t255 * t263) * t281;
t229 = t232 * t297 + t247 * t341;
t228 = t231 * t297 + t246 * t342;
t227 = t230 * t297 + t245 * t343;
t226 = t232 * t341 + t247 * t272;
t225 = t231 * t342 + t246 * t272;
t224 = t230 * t343 + t245 * t272;
t1 = [(-(-t243 * t329 + t252 * t258) * t265 - (-t243 * t332 + t252 * t259) * t262) * t321 + t258 * t322 - t278 * t311 - t258 * t308 + (-(-t241 * t330 + t250 * t256) * t264 - (-t241 * t333 + t250 * t257) * t261) * t324 + t256 * t325 - t277 * t312 - t256 * t309 + (-(-t239 * t331 + t248 * t254) * t263 - (-t239 * t334 + t248 * t255) * t260) * t327 + t254 * t328 - t276 * t313 - t254 * t310 - (-t279 * t299 + t280 * t300) * t338; (-(-t244 * t329 + t253 * t258) * t265 - (-t244 * t332 + t253 * t259) * t262) * t321 + t259 * t322 - t275 * t311 - t259 * t308 + (-(-t242 * t330 + t251 * t256) * t264 - (-t242 * t333 + t251 * t257) * t261) * t324 + t257 * t325 - t274 * t312 - t257 * t309 + (-(-t240 * t331 + t249 * t254) * t263 - (-t240 * t334 + t249 * t255) * t260) * t327 + t255 * t328 - t273 * t313 - t255 * t310 - (t279 * t300 + t280 * t299) * t338; (-(-t226 * t329 + t229 * t258) * t265 - (-t226 * t332 + t229 * t259) * t262) * t321 + t247 * t220 + (-(-t225 * t330 + t228 * t256) * t264 - (-t225 * t333 + t228 * t257) * t261) * t324 + t246 * t219 + (-(-t224 * t331 + t227 * t254) * t263 - (-t224 * t334 + t227 * t255) * t260) * t327 + t245 * t218 + (t223 - t317) * t232 + (t222 - t318) * t231 + (t221 - t319) * t230;];
taucX  = t1;
