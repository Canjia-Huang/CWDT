#
# Copyright (C) 2021 Swift Navigation Inc.
# Contact: Swift Navigation <dev@swift-nav.com>
#
# This source is subject to the license found in the file 'LICENSE' which must
# be be distributed together with this source. All other rights reserved.
#
# THIS CODE AND INFORMATION IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND,
# EITHER EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND/OR FITNESS FOR A PARTICULAR PURPOSE.
#

if(TARGET eigen)
    return()
endif()

add_library(eigen INTERFACE)

if (DEFINED THIRD_PARTY_INCLUDES_AS_SYSTEM AND NOT THIRD_PARTY_INCLUDES_AS_SYSTEM)
    target_include_directories(eigen INTERFACE
            ${PROJECT_SOURCE_DIR}/third_party/eigen/)
else()
    target_include_directories(eigen SYSTEM INTERFACE
            ${PROJECT_SOURCE_DIR}/third_party/eigen/)
endif()