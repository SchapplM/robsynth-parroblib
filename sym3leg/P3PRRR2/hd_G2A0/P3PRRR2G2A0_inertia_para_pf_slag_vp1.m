% Calculate inertia matrix for parallel robot
% P3PRRR2G2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRR2G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:31
% EndTime: 2020-03-09 21:21:31
% DurationCPUTime: 0.37s
% Computational Cost: add. (867->102), mult. (1644->205), div. (405->5), fcn. (963->21), ass. (0->112)
t332 = m(3) * pkin(1);
t309 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t250 = t309 * m(3) + Icges(3,3);
t270 = sin(qJ(3,3));
t276 = cos(qJ(3,3));
t295 = (rSges(3,1) * t276 - rSges(3,2) * t270) * t332;
t238 = t295 + t250;
t284 = 0.1e1 / pkin(2);
t331 = t238 * t284;
t272 = sin(qJ(3,2));
t278 = cos(qJ(3,2));
t294 = (rSges(3,1) * t278 - rSges(3,2) * t272) * t332;
t239 = t294 + t250;
t330 = t239 * t284;
t274 = sin(qJ(3,1));
t280 = cos(qJ(3,1));
t293 = (rSges(3,1) * t280 - rSges(3,2) * t274) * t332;
t240 = t293 + t250;
t329 = t240 * t284;
t254 = sin(qJ(2,3) + qJ(3,3));
t271 = sin(qJ(2,3));
t247 = t271 * pkin(1) + pkin(2) * t254;
t264 = 0.1e1 / t270;
t328 = t247 * t264;
t255 = sin(qJ(2,2) + qJ(3,2));
t273 = sin(qJ(2,2));
t248 = t273 * pkin(1) + pkin(2) * t255;
t265 = 0.1e1 / t272;
t327 = t248 * t265;
t256 = sin(qJ(2,1) + qJ(3,1));
t275 = sin(qJ(2,1));
t249 = t275 * pkin(1) + pkin(2) * t256;
t266 = 0.1e1 / t274;
t326 = t249 * t266;
t325 = t250 * t284;
t324 = t254 * t264;
t323 = t255 * t265;
t322 = t256 * t266;
t267 = legFrame(3,2);
t257 = sin(t267);
t321 = t257 * t264;
t268 = legFrame(2,2);
t258 = sin(t268);
t320 = t258 * t265;
t269 = legFrame(1,2);
t259 = sin(t269);
t319 = t259 * t266;
t260 = cos(t267);
t318 = t260 * t264;
t261 = cos(t268);
t317 = t261 * t265;
t262 = cos(t269);
t316 = t262 * t266;
t285 = 0.1e1 / pkin(1);
t315 = t264 * t285;
t314 = t265 * t285;
t313 = t266 * t285;
t312 = t271 * t270;
t311 = t273 * t272;
t310 = t275 * t274;
t277 = cos(qJ(2,3));
t241 = (pkin(2) * t276 + pkin(1)) * t277 - pkin(2) * t312;
t308 = t241 * t321;
t307 = t241 * t318;
t279 = cos(qJ(2,2));
t242 = (pkin(2) * t278 + pkin(1)) * t279 - pkin(2) * t311;
t306 = t242 * t320;
t305 = t242 * t317;
t281 = cos(qJ(2,1));
t243 = (pkin(2) * t280 + pkin(1)) * t281 - pkin(2) * t310;
t304 = t243 * t319;
t303 = t243 * t316;
t244 = t277 * t276 - t312;
t302 = t244 * t321;
t301 = t244 * t318;
t245 = t279 * t278 - t311;
t300 = t245 * t320;
t299 = t245 * t317;
t246 = t281 * t280 - t310;
t298 = t246 * t319;
t297 = t246 * t316;
t296 = Icges(2,3) + Icges(3,3) + (pkin(1) ^ 2 + t309) * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t235 = 0.2e1 * t295 + t296;
t292 = (t235 * t244 - t241 * t331) * t315;
t236 = 0.2e1 * t294 + t296;
t291 = (t236 * t245 - t242 * t330) * t314;
t237 = 0.2e1 * t293 + t296;
t290 = (t237 * t246 - t243 * t329) * t313;
t289 = (t238 * t244 - t241 * t325) * t315;
t288 = (t239 * t245 - t242 * t325) * t314;
t287 = (t240 * t246 - t243 * t325) * t313;
t263 = m(1) + m(2) + m(3);
t286 = (-t257 * t260 - t258 * t261 - t259 * t262) * t263;
t234 = (-t240 * t256 + t249 * t325) * t313;
t233 = (-t239 * t255 + t248 * t325) * t314;
t232 = (-t238 * t254 + t247 * t325) * t315;
t231 = (-t237 * t256 + t249 * t329) * t313;
t230 = (-t236 * t255 + t248 * t330) * t314;
t229 = (-t235 * t254 + t247 * t331) * t315;
t228 = t262 * t287;
t227 = t261 * t288;
t226 = t260 * t289;
t225 = t259 * t287;
t224 = t258 * t288;
t223 = t257 * t289;
t222 = t262 * t290;
t221 = t261 * t291;
t220 = t260 * t292;
t219 = t259 * t290;
t218 = t258 * t291;
t217 = t257 * t292;
t1 = [m(4) + (t260 ^ 2 + t261 ^ 2 + t262 ^ 2) * t263 + (t217 * t302 + t218 * t300 + t219 * t298 + (-t223 * t308 - t224 * t306 - t225 * t304) * t284) * t285, t286 + (t217 * t301 + t218 * t299 + t219 * t297 + (-t223 * t307 - t224 * t305 - t225 * t303) * t284) * t285, (-t217 * t324 - t218 * t323 - t219 * t322 + (t223 * t328 + t224 * t327 + t225 * t326) * t284) * t285; t286 + (t220 * t302 + t221 * t300 + t222 * t298 + (-t226 * t308 - t227 * t306 - t228 * t304) * t284) * t285, m(4) + (t257 ^ 2 + t258 ^ 2 + t259 ^ 2) * t263 + (t220 * t301 + t221 * t299 + t222 * t297 + (-t226 * t307 - t227 * t305 - t228 * t303) * t284) * t285, (-t220 * t324 - t221 * t323 - t222 * t322 + (t226 * t328 + t227 * t327 + t228 * t326) * t284) * t285; (t229 * t302 + t230 * t300 + t231 * t298 + (-t232 * t308 - t233 * t306 - t234 * t304) * t284) * t285, (t229 * t301 + t230 * t299 + t231 * t297 + (-t232 * t307 - t233 * t305 - t234 * t303) * t284) * t285, m(4) + (-t229 * t324 - t230 * t323 - t231 * t322 + (t232 * t328 + t233 * t327 + t234 * t326) * t284) * t285;];
MX  = t1;
