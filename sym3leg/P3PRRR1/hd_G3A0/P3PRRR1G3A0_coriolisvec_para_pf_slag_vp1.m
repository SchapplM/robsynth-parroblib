% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRRR1G3P3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:07
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRRR1G3P3A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G3P3A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:06:45
% EndTime: 2020-03-09 21:06:46
% DurationCPUTime: 1.26s
% Computational Cost: add. (14082->147), mult. (6774->293), div. (2988->5), fcn. (7908->18), ass. (0->147)
t460 = legFrame(3,2);
t452 = cos(t460);
t457 = pkin(7) + qJ(2,3);
t446 = qJ(3,3) + t457;
t433 = sin(t446);
t436 = cos(t446);
t440 = sin(t457);
t443 = cos(t457);
t532 = 0.1e1 / (t433 * t443 - t440 * t436);
t535 = t452 * t532;
t461 = legFrame(2,2);
t453 = cos(t461);
t458 = pkin(7) + qJ(2,2);
t447 = qJ(3,2) + t458;
t434 = sin(t447);
t437 = cos(t447);
t441 = sin(t458);
t444 = cos(t458);
t531 = 0.1e1 / (t434 * t444 - t441 * t437);
t534 = t453 * t531;
t462 = legFrame(1,2);
t454 = cos(t462);
t459 = pkin(7) + qJ(2,1);
t448 = qJ(3,1) + t459;
t435 = sin(t448);
t438 = cos(t448);
t442 = sin(t459);
t445 = cos(t459);
t530 = 0.1e1 / (t435 * t445 - t442 * t438);
t533 = t454 * t530;
t471 = 0.1e1 / pkin(2);
t515 = t530 * t471;
t516 = t531 * t471;
t517 = t532 * t471;
t501 = t436 * t535;
t465 = xDP(1);
t508 = t465 * t471;
t399 = t501 * t508;
t463 = xDP(3);
t469 = 0.1e1 / pkin(3);
t414 = pkin(2) * t443 + pkin(3) * t436;
t496 = t414 * t535;
t490 = t465 * t496;
t483 = t469 * t490;
t507 = t469 * t414;
t449 = sin(t460);
t464 = xDP(2);
t511 = t449 * t464;
t411 = pkin(2) * t440 + pkin(3) * t433;
t514 = t411 * t469;
t378 = t399 + (-t483 + ((-t433 + t514) * t463 + (-t436 + t507) * t511) * t532) * t471;
t504 = t469 * t471;
t381 = (-t490 + (t411 * t463 + t414 * t511) * t532) * t504;
t529 = t378 * t381;
t499 = t437 * t534;
t400 = t499 * t508;
t415 = pkin(2) * t444 + pkin(3) * t437;
t495 = t415 * t534;
t489 = t465 * t495;
t482 = t469 * t489;
t506 = t469 * t415;
t450 = sin(t461);
t510 = t450 * t464;
t412 = pkin(2) * t441 + pkin(3) * t434;
t513 = t412 * t469;
t379 = t400 + (-t482 + ((-t434 + t513) * t463 + (-t437 + t506) * t510) * t531) * t471;
t382 = (-t489 + (t412 * t463 + t415 * t510) * t531) * t504;
t528 = t379 * t382;
t497 = t438 * t533;
t401 = t497 * t508;
t416 = pkin(2) * t445 + pkin(3) * t438;
t494 = t416 * t533;
t488 = t465 * t494;
t481 = t469 * t488;
t505 = t469 * t416;
t451 = sin(t462);
t509 = t451 * t464;
t413 = pkin(2) * t442 + pkin(3) * t435;
t512 = t413 * t469;
t380 = t401 + (-t481 + ((-t435 + t512) * t463 + (-t438 + t505) * t509) * t530) * t471;
t383 = (-t488 + (t413 * t463 + t416 * t509) * t530) * t504;
t527 = t380 * t383;
t387 = t399 + (-t433 * t463 - t436 * t511) * t517;
t417 = -rSges(3,1) * t443 - t440 * rSges(3,2);
t418 = t440 * rSges(3,1) - rSges(3,2) * t443;
t393 = t417 * t433 + t418 * t436;
t526 = t387 ^ 2 * t393;
t388 = t400 + (-t434 * t463 - t437 * t510) * t516;
t419 = -rSges(3,1) * t444 - t441 * rSges(3,2);
t420 = t441 * rSges(3,1) - rSges(3,2) * t444;
t394 = t419 * t434 + t420 * t437;
t525 = t388 ^ 2 * t394;
t389 = t401 + (-t435 * t463 - t438 * t509) * t515;
t421 = -rSges(3,1) * t445 - t442 * rSges(3,2);
t422 = t442 * rSges(3,1) - rSges(3,2) * t445;
t395 = t421 * t435 + t422 * t438;
t524 = t389 ^ 2 * t395;
t523 = t532 * t411;
t522 = t532 * t449;
t521 = t531 * t412;
t520 = t531 * t450;
t519 = t530 * t413;
t518 = t530 * t451;
t455 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t503 = 0.2e1 * pkin(2) * pkin(3);
t502 = t414 * t522;
t500 = t415 * t520;
t498 = t416 * t518;
t375 = t399 + (-t483 / 0.2e1 + ((-t433 + t514 / 0.2e1) * t463 + (-t436 + t507 / 0.2e1) * t511) * t532) * t471;
t493 = t375 * t381 * t393 * t532;
t376 = t400 + (-t482 / 0.2e1 + ((-t434 + t513 / 0.2e1) * t463 + (-t437 + t506 / 0.2e1) * t510) * t531) * t471;
t492 = t376 * t382 * t394 * t531;
t377 = t401 + (-t481 / 0.2e1 + ((-t435 + t512 / 0.2e1) * t463 + (-t438 + t505 / 0.2e1) * t509) * t530) * t471;
t491 = t377 * t383 * t395 * t530;
t487 = -(rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) - Icges(2,3) - Icges(3,3);
t486 = t436 * t493;
t485 = t437 * t492;
t484 = t438 * t491;
t480 = t433 * t440 + t436 * t443;
t479 = t434 * t441 + t437 * t444;
t478 = t435 * t442 + t438 * t445;
t477 = (t417 * t436 - t418 * t433) * pkin(2);
t476 = (t419 * t437 - t420 * t434) * pkin(2);
t475 = (t421 * t438 - t422 * t435) * pkin(2);
t474 = t480 * pkin(2);
t473 = t479 * pkin(2);
t472 = t478 * pkin(2);
t470 = pkin(2) ^ 2;
t468 = pkin(3) ^ 2;
t439 = t470 + t455;
t423 = -t455 * m(3) - Icges(3,3);
t392 = -Icges(3,3) + (-t455 + t475) * m(3);
t391 = -Icges(3,3) + (-t455 + t476) * m(3);
t390 = -Icges(3,3) + (-t455 + t477) * m(3);
t374 = (pkin(3) * t527 + (t380 * pkin(3) + t389 * t472) * t389) * t515;
t373 = (pkin(3) * t528 + (t379 * pkin(3) + t388 * t473) * t388) * t516;
t372 = (pkin(3) * t529 + (t378 * pkin(3) + t387 * t474) * t387) * t517;
t371 = ((-t478 * t377 * t503 - t380 * t468 - t470 * t389) * t469 * t389 - (pkin(3) + t472) * t527) * t515;
t370 = ((-t479 * t376 * t503 - t379 * t468 - t470 * t388) * t469 * t388 - (pkin(3) + t473) * t528) * t516;
t369 = ((-t480 * t375 * t503 - t378 * t468 - t470 * t387) * t469 * t387 - (pkin(3) + t474) * t529) * t517;
t368 = -t423 * t371 - t392 * t374;
t367 = -t423 * t370 - t391 * t373;
t366 = -t423 * t369 - t390 * t372;
t365 = -t392 * t371 - (t487 + (-t439 + 0.2e1 * t475) * m(3)) * t374;
t364 = -t391 * t370 - (t487 + (-t439 + 0.2e1 * t476) * m(3)) * t373;
t363 = -t390 * t369 - (t487 + (-t439 + 0.2e1 * t477) * m(3)) * t372;
t1 = [(t363 * t501 + t364 * t499 + t365 * t497 + (-t366 * t496 - t367 * t495 - t368 * t494) * t469) * t471 + (0.2e1 * t452 * t486 + 0.2e1 * t453 * t485 + 0.2e1 * t454 * t484 + (t494 * t524 + t495 * t525 + t496 * t526) * t469) * m(3); (-t436 * t363 * t522 - t437 * t364 * t520 - t438 * t365 * t518 + (t366 * t502 + t367 * t500 + t368 * t498) * t469) * t471 + (-0.2e1 * t449 * t486 - 0.2e1 * t450 * t485 - 0.2e1 * t451 * t484 + (-t498 * t524 - t500 * t525 - t502 * t526) * t469) * m(3); (-t433 * t532 * t363 - t434 * t531 * t364 - t435 * t530 * t365 + (t366 * t523 + t367 * t521 + t368 * t519) * t469) * t471 + (-0.2e1 * t433 * t493 - 0.2e1 * t434 * t492 - 0.2e1 * t435 * t491 + (-t519 * t524 - t521 * t525 - t523 * t526) * t469) * m(3);];
taucX  = t1;
