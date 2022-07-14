//
// Created by wzy on 2021/8/10.
//

#include "../../include/DSTree/INodeSegmentSplitPolicy.h"
#include "../../include/DSTree/MeanNodeSegmentSplitPolicy.h"
#include "../../include/DSTree/StdevNodeSegmentSplitPolicy.h"
#include <boost/serialization/export.hpp>

BOOST_CLASS_EXPORT_GUID( MeanNodeSegmentSplitPolicy, "MeanNodeSegmentSplitPolicy" );
BOOST_CLASS_EXPORT_GUID( StdevNodeSegmentSplitPolicy, "StdevNodeSegmentSplitPolicy" );