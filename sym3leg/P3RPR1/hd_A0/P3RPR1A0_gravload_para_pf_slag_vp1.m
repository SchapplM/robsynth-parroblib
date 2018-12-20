% Calculate Gravitation load for parallel robot
% P3RPR1A0
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
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% m [3x1]
%   mass of all robot links (including platform)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 17:54
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taugX = P3RPR1A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:54:27
% EndTime: 2018-12-20 17:54:27
% DurationCPUTime: 0.26s
% Computational Cost: add. (309->86), mult. (492->161), div. (24->3), fcn. (362->14), ass. (0->77)
t478 = legFrame(3,3);
t467 = sin(t478);
t470 = cos(t478);
t455 = -t467 * g(1) + t470 * g(2);
t458 = t470 * g(1) + t467 * g(2);
t475 = rSges(2,3) + qJ(2,3);
t481 = sin(qJ(1,3));
t484 = cos(qJ(1,3));
t487 = pkin(1) + rSges(2,1);
t431 = ((-t455 * t487 - t475 * t458) * m(2) + m(1) * (-rSges(1,1) * t455 + rSges(1,2) * t458)) * t484 + t481 * ((-t455 * t475 + t487 * t458) * m(2) + m(1) * (rSges(1,1) * t458 + rSges(1,2) * t455));
t492 = 0.1e1 / qJ(2,3);
t506 = t431 * t492;
t479 = legFrame(2,3);
t468 = sin(t479);
t471 = cos(t479);
t456 = -t468 * g(1) + t471 * g(2);
t459 = t471 * g(1) + t468 * g(2);
t476 = rSges(2,3) + qJ(2,2);
t482 = sin(qJ(1,2));
t485 = cos(qJ(1,2));
t432 = ((-t456 * t487 - t476 * t459) * m(2) + m(1) * (-rSges(1,1) * t456 + rSges(1,2) * t459)) * t485 + t482 * ((-t456 * t476 + t487 * t459) * m(2) + m(1) * (rSges(1,1) * t459 + rSges(1,2) * t456));
t493 = 0.1e1 / qJ(2,2);
t505 = t432 * t493;
t480 = legFrame(1,3);
t469 = sin(t480);
t472 = cos(t480);
t457 = -t469 * g(1) + t472 * g(2);
t460 = t472 * g(1) + t469 * g(2);
t477 = rSges(2,3) + qJ(2,1);
t483 = sin(qJ(1,1));
t486 = cos(qJ(1,1));
t433 = ((-t457 * t487 - t477 * t460) * m(2) + m(1) * (-rSges(1,1) * t457 + rSges(1,2) * t460)) * t486 + t483 * ((-t457 * t477 + t487 * t460) * m(2) + m(1) * (rSges(1,1) * t460 + rSges(1,2) * t457));
t494 = 0.1e1 / qJ(2,1);
t504 = t433 * t494;
t440 = -t455 * t484 + t458 * t481;
t503 = t440 * t492;
t441 = -t456 * t485 + t459 * t482;
t502 = t441 * t493;
t442 = -t457 * t486 + t460 * t483;
t501 = t442 * t494;
t500 = koppelP(1,1);
t499 = koppelP(2,1);
t498 = koppelP(3,1);
t497 = koppelP(1,2);
t496 = koppelP(2,2);
t495 = koppelP(3,2);
t491 = rSges(3,1);
t490 = rSges(3,2);
t489 = xP(3);
t488 = pkin(1) + pkin(2);
t474 = cos(t489);
t473 = sin(t489);
t466 = t483 * qJ(2,1) + t488 * t486;
t465 = t482 * qJ(2,2) + t488 * t485;
t464 = t481 * qJ(2,3) + t488 * t484;
t463 = -t486 * qJ(2,1) + t483 * t488;
t462 = -t485 * qJ(2,2) + t482 * t488;
t461 = -t484 * qJ(2,3) + t481 * t488;
t454 = -t473 * t497 + t474 * t500;
t453 = -t473 * t496 + t474 * t499;
t452 = -t473 * t495 + t474 * t498;
t451 = -t473 * t500 - t474 * t497;
t450 = -t473 * t499 - t474 * t496;
t449 = -t473 * t498 - t474 * t495;
t448 = -t469 * t483 + t472 * t486;
t447 = t469 * t486 + t472 * t483;
t446 = -t468 * t482 + t471 * t485;
t445 = t468 * t485 + t471 * t482;
t444 = -t467 * t481 + t470 * t484;
t443 = t467 * t484 + t470 * t481;
t439 = -t469 * t463 + t466 * t472;
t438 = -t468 * t462 + t465 * t471;
t437 = -t467 * t461 + t464 * t470;
t436 = t463 * t472 + t469 * t466;
t435 = t462 * t471 + t468 * t465;
t434 = t461 * t470 + t467 * t464;
t1 = [t444 * t506 + t446 * t505 + t448 * t504 - m(3) * g(1) + (-t437 * t503 - t438 * t502 - t439 * t501) * m(2); t443 * t506 + t445 * t505 + t447 * t504 - m(3) * g(2) + (-t434 * t503 - t435 * t502 - t436 * t501) * m(2); m(3) * ((g(1) * t491 + g(2) * t490) * t473 + (g(1) * t490 - g(2) * t491) * t474) + ((t447 * t454 + t448 * t451) * t433 - (t436 * t454 + t439 * t451) * m(2) * t442) * t494 + ((t445 * t453 + t446 * t450) * t432 - (t435 * t453 + t438 * t450) * m(2) * t441) * t493 + ((t443 * t452 + t444 * t449) * t431 - (t434 * t452 + t437 * t449) * m(2) * t440) * t492;];
taugX  = t1;
