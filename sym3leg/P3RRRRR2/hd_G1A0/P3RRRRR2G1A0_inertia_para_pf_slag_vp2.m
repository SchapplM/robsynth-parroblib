% Calculate inertia matrix for parallel robot
% P3RRRRR2G1P1A0
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
%   pkin=[a2,a4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:05
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRRRR2G1P1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR2G1P1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:04:17
% EndTime: 2020-03-09 21:04:18
% DurationCPUTime: 0.67s
% Computational Cost: add. (2001->180), mult. (1926->293), div. (570->14), fcn. (1446->42), ass. (0->167)
t573 = -2 * pkin(1);
t488 = cos(qJ(3,3));
t471 = 0.1e1 / t488;
t572 = t471 ^ 2;
t490 = cos(qJ(3,2));
t473 = 0.1e1 / t490;
t571 = t473 ^ 2;
t492 = cos(qJ(3,1));
t475 = 0.1e1 / t492;
t570 = t475 ^ 2;
t483 = sin(qJ(2,3));
t569 = pkin(1) * t483;
t485 = sin(qJ(2,2));
t568 = pkin(1) * t485;
t487 = sin(qJ(2,1));
t567 = pkin(1) * t487;
t489 = cos(qJ(2,3));
t566 = t489 * pkin(1);
t491 = cos(qJ(2,2));
t565 = t491 * pkin(1);
t493 = cos(qJ(2,1));
t564 = t493 * pkin(1);
t563 = Ifges(3,1) + Ifges(2,3);
t482 = sin(qJ(3,3));
t562 = Ifges(3,4) * t482;
t484 = sin(qJ(3,2));
t561 = Ifges(3,4) * t484;
t486 = sin(qJ(3,1));
t560 = Ifges(3,4) * t486;
t464 = qJ(1,3) + legFrame(3,3);
t452 = qJ(2,3) + t464;
t446 = qJ(3,3) + t452;
t447 = -qJ(3,3) + t452;
t419 = sin(t464) * t573 + (-sin(t447) - sin(t446)) * pkin(2);
t431 = 0.1e1 / (sin(qJ(2,3) + qJ(3,3)) + sin(qJ(2,3) - qJ(3,3)));
t559 = t419 * t431;
t465 = qJ(1,2) + legFrame(2,3);
t453 = qJ(2,2) + t465;
t448 = qJ(3,2) + t453;
t449 = -qJ(3,2) + t453;
t420 = sin(t465) * t573 + (-sin(t449) - sin(t448)) * pkin(2);
t432 = 0.1e1 / (sin(qJ(2,2) + qJ(3,2)) + sin(qJ(2,2) - qJ(3,2)));
t558 = t420 * t432;
t466 = qJ(1,1) + legFrame(1,3);
t454 = qJ(2,1) + t466;
t450 = qJ(3,1) + t454;
t451 = -qJ(3,1) + t454;
t421 = sin(t466) * t573 + (-sin(t451) - sin(t450)) * pkin(2);
t433 = 0.1e1 / (sin(qJ(2,1) + qJ(3,1)) + sin(qJ(2,1) - qJ(3,1)));
t557 = t421 * t433;
t422 = cos(t464) * t573 + (-cos(t447) - cos(t446)) * pkin(2);
t556 = t422 * t431;
t423 = cos(t465) * t573 + (-cos(t449) - cos(t448)) * pkin(2);
t555 = t423 * t432;
t424 = cos(t466) * t573 + (-cos(t451) - cos(t450)) * pkin(2);
t554 = t424 * t433;
t425 = t488 * (-mrSges(3,2) * t569 + Ifges(3,6)) - t482 * (mrSges(3,1) * t569 - Ifges(3,5));
t553 = t425 * t471;
t426 = t490 * (-mrSges(3,2) * t568 + Ifges(3,6)) - t484 * (mrSges(3,1) * t568 - Ifges(3,5));
t552 = t426 * t473;
t427 = t492 * (-mrSges(3,2) * t567 + Ifges(3,6)) - t486 * (mrSges(3,1) * t567 - Ifges(3,5));
t551 = t427 * t475;
t494 = 0.1e1 / pkin(2);
t550 = t431 * t494;
t549 = t432 * t494;
t548 = t433 * t494;
t496 = t488 ^ 2;
t547 = (t488 * pkin(2) + t566) / t496;
t497 = t490 ^ 2;
t546 = (t490 * pkin(2) + t565) / t497;
t498 = t492 ^ 2;
t545 = (t492 * pkin(2) + t564) / t498;
t440 = sin(t452);
t468 = 0.1e1 / t483;
t544 = t440 * t468;
t441 = sin(t453);
t469 = 0.1e1 / t485;
t543 = t441 * t469;
t442 = sin(t454);
t470 = 0.1e1 / t487;
t542 = t442 * t470;
t443 = cos(t452);
t541 = t443 * t468;
t444 = cos(t453);
t540 = t444 * t469;
t445 = cos(t454);
t539 = t445 * t470;
t538 = t468 * t482;
t537 = t469 * t484;
t536 = t470 * t486;
t535 = t471 * t494;
t534 = t473 * t494;
t533 = t475 * t494;
t532 = 0.2e1 * t562;
t531 = 0.2e1 * t561;
t530 = 0.2e1 * t560;
t480 = Ifges(3,2) - Ifges(3,1);
t458 = t480 * t496;
t529 = t458 + t563;
t459 = t480 * t497;
t528 = t459 + t563;
t460 = t480 * t498;
t527 = t460 + t563;
t461 = mrSges(3,1) * t566;
t481 = mrSges(2,2) - mrSges(3,3);
t501 = (-(t482 * mrSges(3,2) - mrSges(2,1)) * t489 - t481 * t483) * pkin(1);
t416 = (t461 + t532) * t488 + t501 + t529;
t526 = t416 * t550;
t462 = mrSges(3,1) * t565;
t500 = (-(t484 * mrSges(3,2) - mrSges(2,1)) * t491 - t481 * t485) * pkin(1);
t417 = (t462 + t531) * t490 + t500 + t528;
t525 = t417 * t549;
t463 = mrSges(3,1) * t564;
t499 = (-(t486 * mrSges(3,2) - mrSges(2,1)) * t493 - t481 * t487) * pkin(1);
t418 = (t463 + t530) * t492 + t499 + t527;
t524 = t418 * t548;
t428 = t488 * t532 + t529;
t523 = t428 * t550;
t429 = t490 * t531 + t528;
t522 = t429 * t549;
t430 = t492 * t530 + t527;
t521 = t430 * t548;
t434 = Ifges(3,5) * t482 + Ifges(3,6) * t488;
t520 = t431 * t434 * t471;
t435 = Ifges(3,5) * t484 + Ifges(3,6) * t490;
t519 = t432 * t435 * t473;
t436 = Ifges(3,5) * t486 + Ifges(3,6) * t492;
t518 = t433 * t436 * t475;
t517 = t482 * t547;
t516 = t494 * t547;
t515 = t484 * t546;
t514 = t494 * t546;
t513 = t486 * t545;
t512 = t494 * t545;
t511 = t471 * t538;
t495 = 1 / pkin(1);
t510 = t495 * t538;
t509 = t473 * t537;
t508 = t495 * t537;
t507 = t475 * t536;
t506 = t495 * t536;
t505 = t434 * t535;
t504 = t435 * t534;
t503 = t436 * t533;
t502 = Ifges(1,3) + ((m(2) + m(3)) * pkin(1) ^ 2) + t563;
t415 = t460 + 0.2e1 * (t463 + t560) * t492 + 0.2e1 * t499 + t502;
t414 = t459 + 0.2e1 * (t462 + t561) * t490 + 0.2e1 * t500 + t502;
t413 = t458 + 0.2e1 * (t461 + t562) * t488 + 0.2e1 * t501 + t502;
t412 = t503 + (t418 * t475 - t430 * t512) * t506;
t411 = t504 + (t417 * t473 - t429 * t514) * t508;
t410 = t505 + (t416 * t471 - t428 * t516) * t510;
t409 = (t418 * t539 + t424 * t521) * t495;
t408 = (t417 * t540 + t423 * t522) * t495;
t407 = (t416 * t541 + t422 * t523) * t495;
t406 = (t418 * t542 + t421 * t521) * t495;
t405 = (t417 * t543 + t420 * t522) * t495;
t404 = (t416 * t544 + t419 * t523) * t495;
t403 = (t415 * t539 + t424 * t524) * t495;
t402 = (t414 * t540 + t423 * t525) * t495;
t401 = (t413 * t541 + t422 * t526) * t495;
t400 = (t415 * t542 + t421 * t524) * t495;
t399 = (t414 * t543 + t420 * t525) * t495;
t398 = (t413 * t544 + t419 * t526) * t495;
t397 = t427 * t533 + (t415 * t475 - t418 * t512) * t506;
t396 = t426 * t534 + (t414 * t473 - t417 * t514) * t508;
t395 = t425 * t535 + (t413 * t471 - t416 * t516) * t510;
t1 = [m(4) + (t401 * t541 + t402 * t540 + t403 * t539 + (t407 * t556 + t408 * t555 + t409 * t554) * t494) * t495, (t401 * t544 + t402 * t543 + t403 * t542 + (t407 * t559 + t408 * t558 + t409 * t557) * t494) * t495, (t401 * t511 + t402 * t509 + t403 * t507 + ((t422 * t520 + t423 * t519 + t424 * t518) * t494 + (-t409 * t513 + t445 * t551) * t470 + (-t408 * t515 + t444 * t552) * t469 + (-t407 * t517 + t443 * t553) * t468) * t494) * t495; (t398 * t541 + t399 * t540 + t400 * t539 + (t404 * t556 + t405 * t555 + t406 * t554) * t494) * t495, m(4) + (t398 * t544 + t399 * t543 + t400 * t542 + (t404 * t559 + t405 * t558 + t406 * t557) * t494) * t495, (t398 * t511 + t399 * t509 + t400 * t507 + ((t419 * t520 + t420 * t519 + t421 * t518) * t494 + (-t406 * t513 + t442 * t551) * t470 + (-t405 * t515 + t441 * t552) * t469 + (-t404 * t517 + t440 * t553) * t468) * t494) * t495; (t395 * t541 + t396 * t540 + t397 * t539 + (t410 * t556 + t411 * t555 + t412 * t554) * t494) * t495, (t395 * t544 + t396 * t543 + t397 * t542 + (t410 * t559 + t411 * t558 + t412 * t557) * t494) * t495, m(4) + (t395 * t511 + t396 * t509 + t397 * t507) * t495 + ((t570 + t571 + t572) * t494 * Ifges(3,3) + ((t427 * t570 + (-t412 - t503) * t545) * t536 + (t426 * t571 + (-t411 - t504) * t546) * t537 + (t425 * t572 + (-t410 - t505) * t547) * t538) * t495) * t494;];
MX  = t1;
